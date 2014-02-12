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
	App (string);

	void run (void);
	void run (string);
};

#endif
