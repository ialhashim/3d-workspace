/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#pragma once

struct kdhyperrect {
	int dim;
	double *min, *max;              /* minimum/maximum coords */
};

struct kdnode {
	double *pos;
	int dir;
	int index;

	struct kdnode *left, *right;	/* negative/positive side */
};

struct res_node {
	struct kdnode *item;
	double dist_sq;
	struct res_node *next;
};

struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};

struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};

#include <float.h>
#define kdEpsilon			DBL_EPSILON

#include <vector>

class KDTree
{

public:
	struct kdtree * tree;

	// Constructor
	KDTree(int k = 3);
	~KDTree();

	/* create a kd-tree for "k"-dimensional data */
	struct kdtree *create(int k);

	/* free the struct kdtree */
	void kd_free();

	/* remove all the elements from the tree */
	void clear();

	/* if called with non-null 2nd argument, the function provided
	 * will be called on data pointers (see kd_insert) when nodes
	 * are to be removed from the tree.
	 */
	void data_destructor( void (*destr)(void*));

	/* insert a node, specifying its position, and optional data */
	int insert( const double *pos, int data);
	int insertf( const float *pos, int data);
	int insert3( double x, double y, double z, int data);
	int insert3f( float x, float y, float z, int data);

	/* Find the nearest node from a given point.
	 *
	 * This function returns a pointer to a result set with at most one element.
	 */
	struct kdres *nearest( const double *pos);
	struct kdres *nearestf( const float *pos);
	struct kdres *nearest3( double x, double y, double z);
	struct kdres *nearest3f( float x, float y, float z);

	struct kdnode *get_nearest(float x, float y, float z);

	bool has( const double *pos, double eps = kdEpsilon );
	bool has( double x, double y, double z, double eps = kdEpsilon );

	struct kdres *get_all();
	std::vector<double*> getAll();

	/* Find the N nearest nodes from a given point.
	 *
	 * This function returns a pointer to a result set, with at most N elements,
	 * which can be manipulated with the kd_res_* functions.
	 * The returned pointer can be null as an indication of an error. Otherwise
	 * a valid result set is always returned which may contain 0 or more elements.
	 * The result set must be deallocated with res_free after use.
	 */
	/*
	struct kdres *nearest_n( const double *pos, int num);
	struct kdres *nearest_nf( const float *pos, int num);
	struct kdres *nearest_n3( double x, double y, double z);
	struct kdres *nearest_n3f( float x, float y, float z);
	*/

	/* Find any nearest nodes from a given point within a range.
	 *
	 * This function returns a pointer to a result set, which can be manipulated
	 * by the kd_res_* functions.
	 * The returned pointer can be null as an indication of an error. Otherwise
	 * a valid result set is always returned which may contain 0 or more elements.
	 * The result set must be deallocated with kd_res_free after use.
	 */
	struct kdres *nearest_range( const double *pos, double range);
	struct kdres *nearest_rangef( const float *pos, float range);
	struct kdres *nearest_range3( double x, double y, double z, double range);
	struct kdres *nearest_range3f( float x, float y, float z, float range);

	/* returns the data pointer (can be null) of the current result set item
	 * and optionally sets its position to the pointers(s) if not null.
	 */
	int res_item(struct kdres *set, double *pos);
	int res_itemf(struct kdres *set, float *pos);
	int res_item3(struct kdres *set, double *x, double *y, double *z);
	int res_item3f(struct kdres *set, float *x, float *y, float *z);

	/* equivalent to res_item(set, 0) */
	int res_item_data(struct kdres *set);
	int getData( const double *pos );
};

/* frees a result set returned by kd_nearest_range() */
void kd_res_free(struct kdres *set);

/* returns the size of the result set (in elements) */
int kres_size(struct kdres *set);

/* rewinds the result set iterator */
void kd_res_rewind(struct kdres *set);

/* returns non-zero if the set iterator reached the end after the last element */
int kd_res_end(struct kdres *set);

/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */
int kd_res_next(struct kdres *set);

void clear_rec(struct kdnode *node, void (*destr)(void*));
int insert_rec(struct kdnode **node, const double *pos, int data, int dir, int dim);
int rlist_insert(struct res_node *list, struct kdnode *item, double dist_sq);
void clear_results(struct kdres *set);

struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max);
void hyperrect_free(struct kdhyperrect *rect);
struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
void hyperrect_extend(struct kdhyperrect *rect, const double *pos);
double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos);
