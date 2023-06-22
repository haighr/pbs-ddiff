void model_parameters::MSYsetup() {
 //****************************************************************************
 // Sets up files and parameters for MSY stuff
 //****************************************************************************

    int Nstrat,ii;
    double i;

    //Nstrat = floor((EndStrategy-StartStrategy)/StepStrategy)+1;
    outMSY << "V0" << "\t" << "SB0";
    ii = 0;
    for(i=StartStrategy;i<=EndStrategy;i=i+StepStrategy) {
        ii++;
        outMSY << "\t" << "nProj_" << ii << "\tYield_" << ii << "\tU_" << ii << "\tVB_" << ii << "\tSB_" << ii;
    }
    outMSY << endl;

    ifstream MSYctl("Yields.ctl");
    if(!MSYctl) {
        cerr << "Unable to open Yields.ctl. Will use defaults." << endl;
        MSYmaxIter = 1000;
        MSYtol = 0.01;
    }
    else {
        while(MSYctl.peek() == '#') {
            while(MSYctl.get() != '\n') {}  //advances past the entire line
        }
        MSYctl >> MSYmaxIter;
        MSYctl.get(); //to advance past end of line
        while(MSYctl.peek() == '#') {
            while(MSYctl.get() != '\n') {}  //advances past the entire line
        }
        MSYctl >> MSYtol;
    }
    MSYctl.close();
}

void model_parameters::MSY()
{
 //********************************************************************************************
 // Calculates MSY
 // outputs "nProj Catch U V0 VB" for each catch level
 //********************************************************************************************

    int nproj,icatch;

    double SBnew, SBold, VBnew, VBold, VB0, diff, theCatch, Nold, Nnew, harvRate;
    dmatrix Nmat(1,Nsexes,1,Nages);

    strategySelect();  //find selectivity, which will be strategyva(sex,age)

    Nmat = determInitCond(VBnew, SBnew); //deterministic initial conditions
    SBold = SBnew;   //same as S[StartYear]
    Nold = Nnew = sum(Nmat);
    //cout << "SB0=  " << value(S[StartYear]) << endl;
    //cout << "SBnew=" << SBnew << endl;

    //femaleN << Nmat[1] << endl;
    //maleN << Nmat[2] << endl;

    theCatch = 0;
    oneYearDetermProj(theCatch,Nmat,VB0,SBnew); //to get VB0 using the strategy selectivity curve   DON'T NEED
    outMSY << VB0 << "\t" << SBnew;
    //cout << "SB0=  " << value(S[StartYear]) << endl;
    //cout << "SBnew=" << SBnew << endl << endl;

    //femaleN << Nmat[1] << endl;
    //maleN << Nmat[2] << endl;

    for(theCatch=StartStrategy;theCatch<=EndStrategy;theCatch=theCatch+StepStrategy) {
        if(StrategyType == 2) {
            harvRate = theCatch; //because working over harvest rates and theCatch is changed below, save for output later
        }
        oneYearDetermProj(theCatch,Nmat,VBnew,SBnew);

        //femaleN << Nmat[1] << endl;
        //maleN << Nmat[2] << endl;

        Nnew = sum(Nmat);
        //diff = fabs(SBold-SBnew);
        diff = fabs(Nold-Nnew);
        nproj=1;
        while(diff > MSYtol & nproj < MSYmaxIter) {
            SBold = SBnew;
            Nold = Nnew;
            if(StrategyType == 2) {
                theCatch = harvRate; //because working over harvest rates and theCatch is changed by reference
            }
            oneYearDetermProj(theCatch,Nmat,VBnew,SBnew);

            //femaleN << Nmat[1] << endl;
            //maleN << Nmat[2] << endl;

            Nnew = sum(Nmat);
            diff = fabs(Nold-Nnew);
            diff = fabs(SBold-SBnew);
            nproj++;
        }
        //now population should be at equilibrium
        //cout << "nproj: " << nproj << "  diff: " << diff << " theCatch: " << theCatch << endl;
        if(StrategyType == 1) {
            if(VBnew>0) {harvRate=theCatch/VBnew;}
            else {harvRate=0;}
        }
        outMSY << "\t" << nproj << "\t" << theCatch << "\t" << harvRate << "\t" << VBnew << "\t" << SBnew;
        Nmat = value(N[StartYear]);
        SBold = SBnew = value(S[StartYear]);
        if(StrategyType == 2) {
            theCatch = harvRate; //because working over harvest rates and theCatch is changed by reference
        }
    }
    outMSY << endl;
}

void model_parameters::strategySelect() {
 //*******************************************************
 //finds the selectivity for projections and yields
 //fills in strategyva, which is a class member
 //dynamics function has to be run to get u
 //*******************************************************
    int method,sex,age;

    strategyva.initialize();
    dvariable sumu=0;
    // selectivity at age for the projections
    for (method=1;method<=Nmethods;method++) {
        for (sex=1;sex<=Nsexes;sex++) {
            for (age=1;age<=Nages;age++) {
                if (proj_gears==0) {
                    strategyva[sex][age]+=strategyuproportion[method]*va[method][EndYear][sex][age];
                }
                else {
                    strategyva[sex][age]+=u[method][EndYear]*va[method][EndYear][sex][age];
                }
            }
        }
        sumu+=u[method][EndYear];
    }
    if (proj_gears!=0) strategyva=strategyva/sumu;
}

void model_parameters::oneYearDetermProj(double& Catch, dmatrix& N2, double& VB, double& SSB)
{
 //*******************************************************
 //projects the population one year for a certain catch or harvest rate
 //  if StrategyType is 2, then it is a harvest rate
 //uses same selectivity as projections (strategySelect function)
 //this function was written by Allan Hicks on 3/DEC/07
 //ARGUMENTS: Catch--the catch to be taken out that year or the harvest rate
 //                  if harvest rate, then it is changed to catch
 //           N2--the sex/age structured population matrix (1:Nsexes,1:Nages)
 //
 //RETURN: N2 matrix by reference
 //        VB (vulnerable biomass) by reference
 //        SSB (spawning biomass, female only) by reference
 //*********************************************************
    int method,sex,age;
    double rec, hr;
    dmatrix wProj(1,Nsexes,1,Nages);
    dvector fecProj(1,Nages);

    // Weight and Fecundity for the projection --> EndYear +1 (from setup_variable_parameters)
    wProj = value(w[EndYear+1]);
    fecProj = value(fec[EndYear+1]);

    rec = value(SSB/(alpha+beta*SSB));
    if(StrategyType == 2) {
        hr = Catch;
        Catch = VB*hr;
    }
    else {
        if(VB>0) { hr = Catch/VB; }
        else { hr=0; }
        if(hr > 1) {hr=1;}
    }
    //cout << "hr: " << hr << endl;
    //if(hr>1) {cout << N2 << endl;}
    SSB = VB = 0;
    for (sex=1;sex<=Nsexes;sex++) {
        //work backwards through ages since overwriting numbers at age
        //plus group
        N2[sex][Nages] = (N2[sex][Nages-1]*(1-hr*value(strategyva[sex][Nages-1]))+N2[sex][Nages]*(1-hr*value(strategyva[sex][Nages])))*value(survVec[sex][Nages]);
        if(N2[sex][Nages]<0) {N2[sex][Nages]=0;}
        VB += N2[sex][Nages]*wProj[sex][Nages]*value(strategyva[sex][Nages])*value(survVec05[sex][Nages]);
        if(sex==1) SSB += N2[1][Nages]*fecProj[Nages];
        for (age=Nages-1;age>=2;age--) {
            N2[sex][age] = N2[sex][age-1]*value(survVec[sex][age-1])*(1-hr*value(strategyva[sex][age-1]));
            if(N2[sex][age]<0) { N2[sex][age]=0; }
            VB += N2[sex][age]*wProj[sex][age]*value(strategyva[sex][age])*value(survVec05[sex][age]);
            if(sex==1) {SSB += N2[1][age]*fecProj[age];}
        }
        //recruitment
        N2[sex][1] = rec*RecFraction[sex];
        if(N2[sex][1]<0) { N2[sex][1]=0; }
        VB += N2[sex][1]*wProj[sex][1]*value(strategyva[sex][1])*value(survVec05[sex][1]);
        if(sex==1) {SSB += N2[1][1]*fecProj[1];}
    } //sex
}
