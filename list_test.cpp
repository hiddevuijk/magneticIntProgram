#include <iostream>
#include "clInfo.h"


using namespace std;

int main()
{

	CLInfo a,b,c;
	a.prev_ptr = 0;
	a.next_ptr = &b;
	b.prev_ptr = &a;
	b.next_ptr = &c;
	c.prev_ptr = &b;
	c.next_ptr = 0;
	a.index = 1;
	b.index = 2;
	c.index = 3;

	++a;
	cout << a.index << endl;




	return 0;
}
