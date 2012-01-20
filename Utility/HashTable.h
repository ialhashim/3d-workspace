#pragma once

#if defined(_HAS_TR1)

	#include <unordered_set>
	#include <unordered_map>

	typedef std::tr1::unordered_set<int> BaseIntSet;
	typedef std::tr1::unordered_set<int>::iterator BaseIntSetIterator;
	typedef std::tr1::unordered_map<int, int> VertexMap;
	#define HashMap std::tr1::unordered_map

#else

	#include <set>
	#include <map>

	typedef std::set<int> BaseIntSet;
	typedef std::set<int>::iterator BaseIntSetIterator;
	typedef std::map<int, int> VertexMap;
	#define HashMap std::map

#endif

struct IntSet
{
	BaseIntSet set;

	IntSet(){	}

	IntSet(const IntSet& fromSet){
		this->set = fromSet.set;
	}

	IntSet& operator= (const IntSet& fromSet){
		this->set = fromSet.set;
		return *this;
	}

	inline void insert(int value){
		set.insert( value );
	}

	inline bool has(int value){
		return set.find(value) != set.end();
	}

	inline void remove(int value){
		set.erase(value);
	}

	inline size_t size(){
		return set.size();
	}

	void addAll(const std::vector<int> & elements)
	{
		for(unsigned int i = 0; i < elements.size(); i++)
			this->insert(elements[i]);
	}

	int first(){
		if(set.size())
			return *set.begin();
		else
			return -1;
	}

	std::vector<int> ToVector(){
		std::vector<int> result(set.size());

		BaseIntSet::iterator it = set.begin();

		int i = 0;

		while(it != set.end()){
			result[i++] = *it;
			it++;
		}

		return result;
	}
};
