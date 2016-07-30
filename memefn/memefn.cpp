#include <memory>
#include <iomanip>
#include <stack>
#include "RNA.h"
#include "ErrorChecker.h"

double ComputeEnergy(const std::string& name, const std::string& seq, const std::string& db) {
  auto strand = std::make_unique<RNA>(seq.c_str());
  auto error = std::make_unique<ErrorChecker<RNA>>(strand.get());
  std::stack<int> s;
  for (int i = 0; i < int(db.size()); ++i) {
    if (db[i] == '(') {
      s.push(i);
    } else if (db[i] == ')') {
      strand->SpecifyPair(s.top() + 1, i + 1);
      s.pop();
    }
  }
  double energy = strand->CalculateFreeEnergy(1, true);
  if (error->isErrorStatus()) exit(1);
  return energy;
}


int main(int argc, char* argv[]) {
  while (1) {
    std::string name, seq, db;
    getline(std::cin, name);
    getline(std::cin, seq);
    getline(std::cin, db);
    if (!std::cin) break;
    std::cout << "Energy: " << std::fixed << std::setprecision(1) << ComputeEnergy(name, seq, db) << std::endl;
  }
}
