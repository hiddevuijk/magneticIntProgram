#ifndef GUARD_CLINFO_H
#define GUARD_CLINFO_H

class CLInfo {
public:
	CLInfo(): prev_ptr(0), next_ptr(0) {}

	CLInfo& last() {
		if(next_ptr == 0) return *this;
		else return next_ptr->last();
	}
	CLInfo& first() {
		if(prev_ptr == 0) return *this;
		else return prev_ptr->first();
	}

	void remove() {
		if( next_ptr != 0 )
			next_ptr->prev_ptr = prev_ptr;
		if( prev_ptr != 0 )
			prev_ptr->next_ptr = next_ptr;
		prev_ptr = 0;
		next_ptr = 0;
	}

	void add_after(CLInfo &element) {
		element.next_ptr = next_ptr;
		element.prev_ptr = this;
		if(next_ptr != 0) 
			next_ptr->prev_ptr = &element;
		next_ptr = &element;
	}

	void add_before(CLInfo &element) {
		element.prev_ptr = prev_ptr;
		element.next_ptr = this;
		if(prev_ptr != 0)
			prev_ptr->next_ptr = &element;
		prev_ptr = &element;
			
	}

	CLInfo *prev_ptr;
	CLInfo *next_ptr;
	unsigned int cell_index;
	unsigned int index;
};









#endif
