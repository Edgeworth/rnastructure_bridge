#include <memory>
#include <iomanip>
#include <stack>
#include "RNA.h"
#include "ParseCommandLine.h"

double ComputeEnergy(const std::string& name, const std::string& seq, const std::string& db) {
  auto strand = std::make_unique<RNA>(seq.c_str());
  std::stack<int> s;
  for (int i = 0; i < int(db.size()); ++i) {
    if (db[i] == '(') {
      s.push(i);
    } else if (db[i] == ')') {
      strand->SpecifyPair(s.top() + 1, i + 1);
      s.pop();
    }
  }
  return strand->CalculateFreeEnergy(1, true);
}


int main(int argc, char* argv[]) {
  auto parser = std::make_unique<ParseCommandLine>("memefn");
  parser->addParameterDescription("dot-bracket list", "Name of file containing a list of dot-brackets");
  parser->parseLine(argc, argv);

  if (parser->isError())
    return 1;

  auto inp = std::ifstream(parser->getParameter(1));
  while (1) {
    std::string name, seq, db;
    getline(inp, name);
    getline(inp, seq);
    getline(inp, db);
    if (!inp) break;
    std::cout << "Energy: " << std::fixed << std::setprecision(1) << ComputeEnergy(name, seq, db) << std::endl;
  }
}

