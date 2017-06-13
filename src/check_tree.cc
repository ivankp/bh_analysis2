// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstring>

#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

int main(int argc, char* argv[]) {
  TFile f(argv[1]);
  if (f.IsZombie()) return 1;

  TTree *tree = nullptr;
  unsigned nentries = 0, nevents = 0;
  Int_t prev_id = -1;

  TIter next(f.GetListOfKeys());
  for (TKey *key; (key=static_cast<TKey*>(next())); ) {
    if (!strcmp(key->GetClassName(),"TTree")) {
      if (tree) { // second tree
        if (!strcmp(tree->GetName(),key->GetName())) continue;
        std::cerr << "Second TTree " << key->GetName() << std::endl;
        return 1;
      }
      tree = static_cast<TTree*>(key->ReadObj());
      if (!tree) { // couldn't read tree
        std::cerr << "Couldn\'t read TTree " << key->GetName() << std::endl;
        return 1;
      }

      nentries = tree->GetEntriesFast();
      if (!nentries) {
        std::cerr << "Zero entries in TTree " << tree->GetName() << std::endl;
        return 1;
      }
      
      TTreeReader reader(tree);
      TTreeReaderValue<Int_t> _id(reader,"id");

      while (reader.Next()) {
        Int_t id = *_id;
        if (id!=prev_id) {
          prev_id = id;
          ++nevents;
        }
      }
      if (!nevents) {
        std::cerr << "Zero events in TTree " << tree->GetName() << std::endl;
        return 1;
      }
    }
  }

  std::cout << tree->GetName() << ' '
            << nentries << ' ' << nevents << std::endl;

  return 0;
}
