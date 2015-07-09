// sort algorithm example
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <stdlib.h>
#include <iomanip>

#include <TROOT.h>
#include <TSystem.h>

//bool myfunction (int i,int j) { return (i<j); }

int PlotTest () {
  int myints[] = {32,71,12,45,26,80,53,33};
  std::vector<int> myvector;
  for(UInt_t i=0;i<8;i++)
    myvector.push_back(myints[i]);
    

  // // using default comparison (operator <):
  sort (myvector.begin(), myvector.end());           //(12 32 45 71)26 80 53 33

  // // using function as comp
  // std::sort (myvector.begin()+4, myvector.end(), myfunction); // 12 32 45 71(26 33 53 80)

  // // using object as comp
  // std::sort (myvector.begin(), myvector.end(), myobject);     //(12 26 32 33 45 53 71 80)

  // print out content:
  std::cout << "myvector contains:";
  for (UInt_t i=0;i<myvector.size();i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  return 0;
}
