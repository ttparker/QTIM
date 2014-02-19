#include <time.h>
#include <fstream>
#include "EffectiveHamiltonian.h"
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
    
	// **************** begin modifiable parameters
	const int numberOfTrials = 1;
	int lSys = 16;							// system length - must be even
	std::vector<double> couplingConsts = {-1.};
    double h = 1.;
	int mMax = 16,						    // max number of stored states
		nSweeps = 2;						// number of sweeps to be performed
    double groundStateErrorTolerance = 1e-6;
    
	bool energyOnly = false,
		// calculate observables? (if true, rest of this section irrelevant)
		 calcOneSiteExpValues = true, // calculate single-site expectation values?
		 calcTwoSiteExpValues = true; // calculate two-site expectation values?
	int rangeOfObservables;	// number of sites in center at which to
        // calculate observables (must be even and less than lSys) - set below
	MatrixDd oneSiteOp,
			 firstTwoSiteOp,
			 secondTwoSiteOp;
	if(!energyOnly)
	{
		rangeOfObservables = 8;
        #include "ObservableOps.h"
		oneSiteOp = sigmaz,
		firstTwoSiteOp = sigmax,
		secondTwoSiteOp = sigmax;
	};
	// **************** end modifiable parameters
    
    filein.close();
    
	Hamiltonian ham;			// initialize the system's Hamiltonian
	std::ofstream fileout;
	fileout.open("Output/Output", std::ios::out);
	if (!fileout)
	{
		std::cout << "Couldn't open output file." << std::endl;
		exit(1);
	};
	for(int trial = 0; trial < numberOfTrials; trial++)
	{
		std::cout << "Trial " << trial << ":" <<std::endl;
		fileout << "Trial " << trial << ":" <<std::endl;
		ham.setParams(lSys, couplingConsts, h);
		int lSFinal = ham.lSys / 2 - 1;		// final length of the system block
		std::vector<TheBlock> blocks(ham.lSys - 3);		// initialize system
		blocks[0] = TheBlock(ham, mMax);	// initialize the one-site block
        int skips = 0;
        for(int runningKeptStates = d * d; runningKeptStates <= mMax; skips++)
            runningKeptStates *= d; // find how many edge sites can be skipped
        std::cout << "Performing iDMRG..." << std::endl;
        for(int site = 0; site < skips; site++)
            blocks[site + 1] = blocks[site].nextBlock(ham);       // initial ED
        TheBlock::lancTolerance = groundStateErrorTolerance
                                  * groundStateErrorTolerance / 2;
        for(int site = skips, end = lSFinal - 1; site < end; site++)
            blocks[site + 1] = blocks[site].nextBlock(ham, false);
        if(nSweeps != 0)
            std::cout << "Performing fDMRG..." << std::endl;
        for(int i = 1; i <= nSweeps; i++)			// perform the fDMRG sweeps
		{
            for(int site = lSFinal - 1, end = ham.lSys - 4 - skips; site < end; site++)
                blocks[site + 1] = blocks[site].nextBlock(ham, false, false,
                                                          blocks[ham.lSys - 4 - site],
                                                          blocks[ham.lSys - 5 - site]);
            blocks[skips].reflectPredictedPsi();
                               // reflect the system to reverse sweep direction
            for(int site = skips, end = lSFinal - 1; site < end; site++)
                blocks[site + 1] = blocks[site].nextBlock(ham, false, false,
                                                          blocks[ham.lSys - 4 - site],
                                                          blocks[ham.lSys - 5 - site]);
            std::cout << "Sweep " << i << " complete." << std::endl;
        };
        EffectiveHamiltonian hSuperFinal = blocks[lSFinal - 1]
                                 .createHSuperFinal(ham, skips);
                                               // calculate ground-state energy
		fileout << "Ground state energy density = "
				<< hSuperFinal.gsEnergy / ham.lSys << std::endl	<< std::endl;
		if(!energyOnly)
		{
			std::cout << "Calculating observables..." << std::endl;
			if(calcOneSiteExpValues)   // calculate one-site expectation values
				oneSiteExpValues(oneSiteOp, rangeOfObservables, ham.lSys,
								 hSuperFinal, blocks, fileout);
			if(calcTwoSiteExpValues)   // calculate two-site expectation values
				twoSiteExpValues(firstTwoSiteOp, secondTwoSiteOp,
								 rangeOfObservables, ham.lSys, hSuperFinal,
								 blocks, fileout);
		};
		std::cout << std::endl;
		fileout << std::endl;
	};

	clock_t stop = clock();
    std::cout << "Done." << std::endl;
	fileout << "Elapsed time (s): " << (double)(stop - start)/CLOCKS_PER_SEC
			<< std::endl;
	fileout.close();

	return 0;
}
