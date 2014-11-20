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
*  @brief  Gibt das uebergebene Array aus.
*  @param  __arr  Array
*  @param  __n    Groeße des Arrays
*/
void print(bool* arr, size_t n) {
	std::cout << "\n[ ";
	for (size_t i = 0; i < n; ++i) {
		std::cout << arr[i] << " ";
	}
	std::cout << "]";
}

/**
*  @brief  Gibt das naechstfreie false-Bit zurueck.
*  @param  __arr    Array
*  @param  __n      Groeße des Arrays
*  @param  __index  Startindex (optional)
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
		for (size_t i = 0; i < realEnd; ++i) {
			size_t j = i + lvl + 1;
			if (!(tree & (1 << i)) && !(tree & (1 << j)))
			{
				size_t treeTree = (1 << i) + (1 << j);
				LangfordRecursiveBit( tree + treeTree, lvl - 1, size, count);
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
		// Alle möglichen Positionen für Index i und i + Abstand testen
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
*/
tbb::task* Langford::execute() {
	if (lvl < 1) {
		++count;
	}
	else {
		tbb::task_list taskList;
		size_t taskCount = 0;
		size_t realEnd = size - lvl - 1;
		for (size_t i = 0; i < realEnd; ++i)
		{
			size_t j = (i + lvl + 1);
			if (!(tree & (1 << i)) && !(tree & (1 << j)))
			{
				size_t treeTree = (1 << i) + (1 << j);
				taskList.push_back(*new (allocate_child()) Langford(tree + treeTree, lvl - 1, size, count));
				++taskCount;
			}
		}
		set_ref_count(++taskCount);
		spawn_and_wait_for_all(taskList);
	}
	return NULL;
}
