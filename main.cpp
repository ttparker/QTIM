#include <time.h>
#include <fstream>
#include "FreeFunctions.h"

using namespace Eigen;

int main()
{
    clock_t start = clock();
    
    std::ifstream filein("Input/Input");
    if(!filein)
    {
        std::cerr << "Couldn't open input file." << std::endl;
        exit(EXIT_FAILURE);
    };
    
	// read in parameters that are constant across all trials:
    int nTrials,
        nCouplingConstants;
    bool calcObservables;
	filein >> nTrials >> nCouplingConstants >> calcObservables;
    bool calcOneSiteExpValues;
    int indexOfOneSiteOp;
    MatrixDd oneSiteOp;
    bool calcTwoSiteExpValues;
    int indexOfFirstTwoSiteOp,
        indexOfSecondTwoSiteOp;
    MatrixDd firstTwoSiteOp,
             secondTwoSiteOp;
    #include "ObservableOps.h"
    if(calcObservables)
    {
        filein >> calcOneSiteExpValues;
        if(calcOneSiteExpValues)
        {
            filein >> indexOfOneSiteOp;
            oneSiteOp = obsList[indexOfOneSiteOp];
        };
        filein >> calcTwoSiteExpValues;
        if(calcTwoSiteExpValues)
        {
            filein >> indexOfFirstTwoSiteOp >> indexOfSecondTwoSiteOp;
            firstTwoSiteOp = obsList[indexOfFirstTwoSiteOp];
            secondTwoSiteOp = obsList[indexOfSecondTwoSiteOp];
        };
    };
    
    std::ofstream infoout("Output/Info");
    if(infoout)
        infoout << "Number of trials: " << nTrials
                << "\nCalculate one-site observables? "
                << (calcOneSiteExpValues ? "Yes" : "No")
                << "\nIndex of one-site observable: " << indexOfOneSiteOp
                << "\nCalculate two-site observables? "
                << (calcTwoSiteExpValues ? "Yes" : "No")
                << "\nIndices of two-site observables: " << indexOfFirstTwoSiteOp
                << " " << indexOfSecondTwoSiteOp << std::endl;
    else
    {
        std::cerr << "Couldn't open output files." << std::endl;
        exit(EXIT_FAILURE);
    };
    Hamiltonian ham;			// initialize the system's Hamiltonian
	for(int trial = 1; trial <= nTrials; trial++)
	{
        clock_t startTrial = clock();
		std::cout << "Trial " << trial << ":" <<std::endl;
        std::ofstream fileout("Output/Trial_" + std::to_string(trial));
		fileout << "Trial " << trial << ":\n" << std::endl;
        
        // read in parameters that vary over trials:
        int lSys;                           // system length
        filein >> lSys;
        std::vector<double> couplingConstants(nCouplingConstants);
        for(int i = 0; i < nCouplingConstants; i++)
            filein >> couplingConstants[i];
        int rangeOfObservables;
        double groundStateErrorTolerance;
        int mMax,                           // max number of stored states
            nSweeps;                        // number of sweeps to be performed
        filein >> rangeOfObservables >> groundStateErrorTolerance >> mMax
               >> nSweeps;
        if(rangeOfObservables == -1)
            rangeOfObservables = lSys;
        
        fileout << "System length: " << lSys << "\nCoupling constants:";
        for(int i = 0; i < nCouplingConstants; i++)
            fileout << " " << couplingConstants[i];
        fileout << "\nLanczos tolerance: " << groundStateErrorTolerance
                << "\nBond dimension: " << mMax << "\nNumber of sweeps: "
                << nSweeps << std::endl << std::endl;
        ham.setParams(lSys, couplingConstants);
        int skips = 0;
        for(int runningKeptStates = d * d; runningKeptStates <= mMax; skips++)
            runningKeptStates *= d;  // find how many edge sites can be skipped
        std::vector<TheBlock> leftBlocks(lSys - 3 - skips), // initialize system
                              rightBlocks(lSys - 3 - skips);
		leftBlocks[0] = rightBlocks[0]
                      = TheBlock(ham, mMax);	   // initialize the one-site block
        std::cout << "Performing iDMRG..." << std::endl;
        for(int site = 0; site < skips; site++)
            rightBlocks[site + 1] = leftBlocks[site + 1] 
                                  = leftBlocks[site].nextBlock(ham); // initial ED
        TheBlock::lancTolerance = groundStateErrorTolerance
                                  * groundStateErrorTolerance / 2;
        int lSFinal = lSys / 2 - 1;         // final length of the system block
        for(int site = skips, end = lSFinal - 1; site < end; site++)
        {
            rightBlocks[site + 1] = leftBlocks[site + 1]
                                  = leftBlocks[site].nextBlock(ham, false);
            rightBlocks[site].primeToRhoBasis = leftBlocks[site].primeToRhoBasis;
        };
        if(nSweeps == 0)
            leftBlocks[lSFinal - 1].randomSeed();
        else
        {
            std::cout << "Performing fDMRG..." << std::endl;
            for(int i = 1; i <= nSweeps; i++)			// perform the fDMRG sweeps
            {
                for(int site = lSFinal - 1, end = lSys - 4 - skips; site < end;
                    site++)
                    leftBlocks[site + 1] = leftBlocks[site].nextBlock(ham,
                                           false, false,
                                           rightBlocks[lSys - 4 - site],
                                           rightBlocks[lSys - 5 - site]);
                rightBlocks[skips].reflectPredictedPsi();
                               // reflect the system to reverse sweep direction
                for(int site = skips, end = lSys - 4 - skips; site < end; site++)
                    rightBlocks[site + 1] = rightBlocks[site].nextBlock(ham,
                                            false, false,
                                            leftBlocks[lSys - 4 - site],
                                            leftBlocks[lSys - 5 - site]);
                leftBlocks[skips].reflectPredictedPsi();
                for(int site = skips, end = lSFinal - 1; site < end; site++)
                    leftBlocks[site + 1] = leftBlocks[site].nextBlock(ham,
                                           false, false,
                                           rightBlocks[lSys - 4 - site],
                                           rightBlocks[lSys - 5 - site]);
                std::cout << "Sweep " << i << " complete." << std::endl;
            };
        };
        EffectiveHamiltonian hSuperFinal = leftBlocks[lSFinal - 1]
                      .createHSuperFinal(ham, rightBlocks[lSFinal - 1], skips);
                                               // calculate ground-state energy
		fileout << "Ground state energy density = "
				<< hSuperFinal.gsEnergy / lSys << std::endl	<< std::endl;
		if(calcObservables)
		{
			std::cout << "Calculating observables..." << std::endl;
			if(calcOneSiteExpValues)   // calculate one-site expectation values
				oneSiteExpValues(oneSiteOp, rangeOfObservables, lSys,
								 hSuperFinal, leftBlocks, rightBlocks, fileout);
			if(calcTwoSiteExpValues)   // calculate two-site expectation values
				twoSiteExpValues(firstTwoSiteOp, secondTwoSiteOp,
								 rangeOfObservables, lSys, hSuperFinal,
                                 leftBlocks, rightBlocks, fileout);
		};
        std::cout << std::endl;
        clock_t stopTrial = clock();
		fileout << "Elapsed time: "
                << double(stopTrial - startTrial)/CLOCKS_PER_SEC << " s"
                << std::endl;
        fileout.close();
	};
    filein.close();
    
    clock_t stop = clock();
    std::cout << "Done. Elapsed time: " << double(stop - start)/CLOCKS_PER_SEC
			  << " s" << std::endl;

	return 0;
}
