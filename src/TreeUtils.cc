#include "interface/TreeUtils.h"

std::map<std::string,TLeaf*> initTreeVars(TChain* chain)
{
	std::map<std::string,TLeaf*> map;

	TObjArray* array = chain -> GetListOfBranches();
	int nBranches = array -> GetEntriesFast();

	for(int i=0; i<nBranches; ++i)
	{
		TString branchName = array->At(i)->GetName();

		TLeaf* leaf = chain->GetLeaf(branchName);
		map[(std::string)branchName] = leaf;
	}
	return map;
}

