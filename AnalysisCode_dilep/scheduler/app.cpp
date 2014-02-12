#include "app.h"
#include <iostream>

App::App (string _name, string _args) {
	name = _name;
	args = _args;
}

App::App (string _name) {
	name = _name;
}

void App::run (void) {
	#ifdef VERBOSE
	string command = string("time ./") + name + string(" ") + args;
	#else
	string command = string("./") + name + string(" ") + args;
	#endif
	
	system (command.c_str());
}

void App::run (string _arg) {
	#ifdef VERBOSE
	string command = string("time ./") + name + string(" ") + _arg;
	#else
	string command = string("./") + name + string(" ") + _arg;
	#endif
	
	system (command.c_str());
}
