#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"

#define nLep 8
#define nJet 10



 
std::map<std::string,TLeaf*> initTreeVars(TChain* chain);


#endif
