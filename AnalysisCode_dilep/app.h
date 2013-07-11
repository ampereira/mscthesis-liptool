#ifndef APP
#define APP

#include <cstdlib>
#include <string>

using namespace std;

class App {

	string name;
	string args;

public:
	App (string, string);

	void run (void);

};

#endif
