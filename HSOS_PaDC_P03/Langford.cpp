//============================================================================
// Name        : Langford.cpp
// Author      : Tobias Sibera <Tobias.Sibera@HS-Osnabrueck.de>
// Version     : 1.00
// Copyright   : GPLv3
// Created on  : 17.11.2014
// Description : Algorithmus: Langford. Parallelisierung mit tbb.
//============================================================================

#include "Langford.h"

/**
*  @brief  Gibt ein unsigned long-Wert als binäre Folge aus.
*  @param  x     Array
*  @param  size  Groeße des Arrays
*/
void printBits(size_t& x, size_t& size) {
	size_t half = size >> 1;
    for (size_t i = 0; i < size; ++i) {
    	if (i == half) {
    		std::cout << ' ';
    	}
        std::cout << (x & ((size_t) 1 << (half - i)) ? '1' : '0');
    }
    std::cout << " (" << x << ')';
}

/**
*  @brief  Gibt ein umgekehrte, binären Wert zurück.
*  @param  x     Wert
*  @param  size  Groeße des Arrays (optional)
*  @return Umgekehrte Bitfolge.
*/
size_t reverse(size_t x, size_t size = 0) {
    x = ((x >>  1) & 0x55555555u) | ((x & 0x55555555u) <<  1);
    x = ((x >>  2) & 0x33333333u) | ((x & 0x33333333u) <<  2);
    x = ((x >>  4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) <<  4);
    x = ((x >>  8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) <<  8);
    x = ((x >> 16) & 0xffffu)	  | ((x & 0xffffu) 	   << 16);
    return x >> (32 - size);
}

/**
*  @brief  Gibt das naechstfreie false-Bit zurueck.
*  @param  __arr    Array
*  @param  __n      Groeße des Arrays
*  @param  __index  Startindex (optional)
*  @return Position des nächst falschen Bits.
*/
size_t nextClearBit(bool* arr, size_t n, size_t index = 0) {
	for (size_t i = index; i < n; ++i) {
		if (!arr[i]) {
			return i;
		}
	}
	return n;
}

/**
*  @brief  Langford-Algorithmus (Bit-Shiftig).
*  @param  __arr    Boolsches Array
*  @param  __lvl    SB-Ebene
*  @param  __size   Groeße
*  @param  __count  Globaler zaehler
*/
void LangfordRecursiveBit(size_t tree, size_t lvl, size_t& size, size_t& count) {
	if (lvl < 1) {
		++count;
	}
	else {
		size_t realEnd = size - lvl - 1;
		std::vector<size_t> candidates(0);
		for (size_t i = 0; i < realEnd; ++i) {
			size_t j = i + lvl + 1;
			if (!(tree & (1 << i)) && !(tree & (1 << j))) {
				bool symmetric = false;
				size_t newTree = tree + (1 << i) + (1 << j);
				size_t revTree = reverse(newTree, size);
				for (size_t k = 0; k < candidates.size(); ++k)  {
					if (candidates[k] == revTree) {
						symmetric = true;
						break;
					}
				}
				if (!symmetric) {
					candidates.push_back(newTree);
					LangfordRecursiveBit(newTree, lvl - 1, size, count);
				}
			}
		}
	}
}

/**
*  @brief  Langford-Algorithmus.
*  @param  __arr    Boolsches Array
*  @param  __lvl    SB-Ebene
*  @param  __size   Groeße
*  @param  __count  Globaler zaehler
*/
void LangfordRecursive(bool* arr, size_t lvl, size_t& size, size_t& count) {
	if (lvl < 1) {
		for (size_t i = 0; i < size; ++i) {
			if (!arr[i]) {
				return;
			}
		}
		++count;
	}
	else {
		for (size_t i = -1; (i = nextClearBit(arr, size, i + 1)) < size - lvl - 1;) {
			size_t j = i + lvl + 1;
			if (!arr[j]) {
				arr[i] = !arr[i];
				arr[j] = !arr[j];
				LangfordRecursive(arr, lvl - 1, size, count);
				arr[i] = !arr[i];
				arr[j] = !arr[j];
			}
		}
	}
}

/**
*  @brief  Langford-Algorithmus. Lösung mithilfe von Tasks.
*  @param  __tree   Bit-Sequenz
*  @param  __lvl    SB-Ebene
*  @param  __size   Groeße
*  @param  __count  Globaler zaehler
*/
Langford::Langford(size_t __tree, size_t __lvl, size_t& __size, tbb::atomic<size_t>& __count) : tree(__tree), lvl(__lvl), size(__size), count(__count) {
}

/**
*  @brief  Vererbte und ueberschriebene Methode der tbb::task-Klasse,
*  welche die erbende Klasse in einem Task ausfuehren laesst.
*  @return tbb::task.
*/
tbb::task* Langford::execute() {
	if (lvl < 1) {
		++count;
	}
	else {
		tbb::task_list taskList;
		size_t taskCount = 0;
		size_t realEnd = size - lvl - 1;
		std::vector<size_t> candidates(realEnd >> 1);
		for (size_t i = 0; i < realEnd; ++i) {
			size_t j = i + lvl + 1;
			if (!(tree & (1 << i)) && !(tree & (1 << j))) {
				bool symmetric = false;
				size_t newTree = tree + (1 << i) + (1 << j);
				size_t revTree = reverse(newTree, size);
				candidates.push_back(newTree);
				for (size_t j = 0; j < candidates.size() - 1; ++j)  {
					if (candidates[j] == revTree) {
						symmetric = true;
						break;
					}
				}
				if (!symmetric) {
					candidates.push_back(newTree);
					taskList.push_back(*new (allocate_child()) Langford(newTree, lvl - 1, size, count));
					++taskCount;
				}
			}
		}
		candidates.clear();
		candidates.shrink_to_fit();
		if (taskCount > 0) {
			set_ref_count(++taskCount);
			spawn_and_wait_for_all(taskList);
		}
	}
	return NULL;
}
