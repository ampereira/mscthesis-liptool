#include "app.h"
#include <iostream>

App::App (string _name, string _args) {
	name = _name;
	args = _args;
}

void App::run (void) {
	string command = string("time ./") + name + string(" ") + args;
	//cout << command << endl;
	system (command);
}
