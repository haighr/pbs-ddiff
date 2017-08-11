adstring stripExtensionNew(adstring fileName){
  /*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
  */
  const int length = fileName.size();
  for (int i=length; i>=0; --i){
    if (fileName(i)=='.'){
      return fileName(1,i-1);
    }
  }
  return(fileName);
}

ivector getIndex(const dvector& a, const dvector& b){
  /*
    Returns a dvector which contains the elements of dvector a where
    the elements of b are not NA
    Assumes the length of a and b are the same.
   */
  int i,j,n;
  n = 0;
  j = 1;
  for(i = a.indexmin(); i <= a.indexmax(); i++){
    if(b(i) != NA){
      n++;
    }
  }
  ivector tmp(1,n);
  for(i = a.indexmin(); i <= a.indexmax(); i++){
    if(b(i) != NA){
      tmp(j++) = a(i);
    }
  }
  return(tmp);
}

void write_proj_headers(ofstream &ofsP, int syr, int nyr, int pyr){
  // Write the decision table headers for projection years
  int i;
  ofsP<<"TAC"                   <<",";
  for (i=nyr+1; i<=pyr; i++){
//	for (i=nyr-5; i<=pyr; i++){
    ofsP<<"B"<<i                <<",";
  }
  ofsP<<"B0"                    <<",";
  ofsP<<"B040"                  <<",";
  ofsP<<"B020"                  <<",";
  ofsP<<"B"<<syr                <<",";
  ofsP<<"B"<<nyr+2<<"B"<<nyr+1  <<",";  //want probability B2016<B2015 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"B040"      <<",";  //want probability B2016<0.4B0 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"B020"      <<",";  //want probability B2016<0.4B0 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"B"<<syr    <<",";
  for (i=nyr; i<=(pyr-1); i++){
    ofsP<<"F"<<i                <<",";
  }
  ofsP<<"F"<<nyr+1<<"F"<<nyr    <<",";  //want probability F2015>F2014     - this will be > 1 if true
  ofsP<<"U"<<nyr+1              <<",";
  ofsP<<"U"<<nyr+1<<"U"<<nyr    <<",";
  //MSY based ref points
  ofsP<<"MSY"                   <<",";
  ofsP<<"BMSY"                  <<",";
  ofsP<<"B"<<nyr+2<<"BMSY"      <<",";  //want probability B2016<BMSY - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"BMSY80"    <<",";  //want probability B2016<0.8BMSY - this will be< 1 if true
  ofsP<<"B"<<nyr+2<<"BMSY40"    <<",";  //want probability B2016<0.4BMSY - this will be < 1 if true
  ofsP<<"FMSY"                  <<",";
  ofsP<<"F"<<nyr+1<<"FMSY"      <<",";
  ofsP<<"UMSY"                  <<",";
  ofsP<<"U"<<nyr+1<<"UMSY"      <<endl; //want probability F2015>FMSY - this will be > 1 if true
}

void write_proj_output(ofstream &ofsP, int syr, int nyr, double tac, int pyr, dvector p_sbt, dvector p_ft, double bo, double fmsy, double bmsy, double msy){
	// Write the projection output to the file
	int i;
	ofsP<<tac                    <<",";
	for (i=nyr+1; i<=pyr; i++){
//	for (i=nyr-5; i<=pyr; i++){
		ofsP<<p_sbt(i)            <<",";
	}
//      <<p_sbt(pyr)            <<","
//      <<p_sbt(pyr+1)          <<","
	ofsP<<bo                     <<","
	<<0.4*bo                     <<","
	<<0.2*bo                     <<","
	<<p_sbt(syr)                 <<","
	<<p_sbt(nyr+2)/p_sbt(nyr+1)  <<","
	<<p_sbt(nyr+2)/(0.4*bo)      <<","
	<<p_sbt(nyr+2)/(0.2*bo)      <<","
	<<p_sbt(nyr+2)/p_sbt(syr)    <<",";
	for (i=nyr; i<=(pyr-1); i++){
		ofsP<<p_ft(i)             <<",";
	}
	ofsP<<p_ft(nyr+1)/p_ft(nyr)  <<","
	<<(1. - mfexp(-p_ft(nyr+1))) <<","
	<<(1. - mfexp(-p_ft(nyr+1)))/(1. - mfexp(-p_ft(nyr))) <<","
	//MSY based ref points
	<<msy                        <<","
	<<bmsy                       <<","
	<<p_sbt(nyr+2)/bmsy          <<","
	<<p_sbt(nyr+2)/(0.8*bmsy)    <<","
	<<p_sbt(nyr+2)/(0.4*bmsy)    <<","
	<<fmsy                       <<","
	<<p_ft(nyr+1)/fmsy           <<","
	<<(1. - mfexp(-fmsy))        <<","
	<<(1. - mfexp(-p_ft(nyr+1)))/(1. - mfexp(-fmsy)) <<endl;
}

