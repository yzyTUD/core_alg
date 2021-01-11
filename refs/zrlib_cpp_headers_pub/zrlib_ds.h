#pragma once

#include<iostream>
#include<queue>
#include<stack>
#include<cstring>
#include<string>
#include<complex> // for using complex numbers
#include<cfloat>
#include <random>
#include <unordered_map>
#include <vector>

//#include "zrlib_util.h" // include util package in my lib.

template <typename T> struct sample_set {
	typedef int element_index;
	std::vector<T> elements;
	std::unordered_map<T, element_index> index_lut;

	// create empty set
	sample_set() {}

	// reserve memory for n elements
	void reserve(size_t n) {
		elements.reserve(n);
		index_lut.rehash(n);
	}

	// insert element elem into set
	void insert(const T& elem) {
		// guard against duplicates
		if (index_lut.find(elem) == index_lut.end()) {
			element_index idx = (element_index)elements.size();
			elements.push_back(elem);
			index_lut[elem] = idx;
		}
	}

	// remove element elem from set
	bool remove(const T& elem) {
		auto it = index_lut.find(elem);
		if (it == index_lut.end()) {
			return false;
		}

		int i = it->second;

		std::swap(elements[i], elements.back());
		elements.pop_back();

		index_lut.erase(it);
		if (unsigned(i) < elements.size()) {
			index_lut[elements[i]] = i;
		}
		return true;
	}

	// draw a sample from set
	template <typename Engine> const T& sample(Engine& eng) {
		int b = (element_index)(elements.size() - 1);
		std::uniform_int_distribution<int> uniform_dist(0, b);
		element_index idx = uniform_dist(eng);
		return elements[idx];
	}

	// returns number of elements in set
	size_t size() const { return elements.size(); }

	// returns true if set is empty
	bool empty() const { return elements.empty(); }
};