#include "app.h"
#include <iostream>

App::App (string _name, string _args) {
	name = _name;
	args = _args;
}

void App::run (void) {
	#ifdef VERBOSE
	string command = string("time ./") + name + string(" ") + args;
	#else
	string command = string("./") + name + string(" ") + args;
	#endif
	
	system (command.c_str());
}
