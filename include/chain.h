void chain(Int_t nEvents = 1000000000000, Char_t *inputFile="pico.list", Char_t *outputFile="test.root")
{

  //read the list of files                                                                                                                                                         
  TList nameList;
  char str[2000];
  fstream file_op(inputFile,ios::in);
  while(file_op >> str){
    nameList.Add(new TObjString(str));
  }

  //put the list into chain                                                                                                                                            
  TChain *t = new TChain("PicoDst");
  TIter next(&nameList);
  TObjString *fileNm;
  int i=0;
  while(fileNm=( TObjString*)next()){
    t->AddFile(fileNm->String());
    cout<<i<<": "<<fileNm->String()<<endl;
    i++;
  }

  //load compiled macro 
  gSystem->Load("PicoDst_C.so");
  //create an instance of the PicoDst object
  PicoDst *T = new PicoDst(t,outputFile);
  //run Loop method (see PicoDst.C and PicoDst.h)
  T->Loop();
}
