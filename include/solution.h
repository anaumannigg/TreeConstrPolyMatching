#ifndef _solution_included_
#define _solution_included_

#include "cgal_includes.h"

//class to store a many to many matching as match indices per polygon of both sets
//unmatched polygons get a negative matching index
class Solution {
	std::vector < std::pair<std::vector<int>, std::vector<int>>> matching;
	//assigns every osm polygon an index of the group it belongs to
	std::vector<int> set1_match_index;
	//assigns every osm polygon a match weight corresponding to the weight of the selected edge in E3 by the ILP
	std::vector<double> set1_match_weight;
	//assigns every atkis polygon an index of the group it belongs to
	std::vector<int> set2_match_index;
	//assigns every atkis polygon a match weight corresponding to the weight of the selected edge in E3 by the ILP
	std::vector<double> set2_match_weight;
	//remember amount of matches
	int match_count;
	//remember current match ids (minimum and maximum)
	std::pair<int, int> match_ids;
	//remember target value of solution
	double target_value;
	//remember distribution of matches;
	int numUnmatched, num1to1, num1toN, numMtoN;
	double obj1to1, obj1toN, objMtoN;

public:
    //standard constructor, initializes empty solution
    Solution();

	//init solution instance
	Solution(int num_polys1,int num_polys2);

	//init solution instance with specific beginning IDs for matchings, to work on subsets and merge solution afterwards
	Solution(int num_polys1, int num_polys2, std::pair<int, int> match_id_limits);

	void addMatch(std::vector<int> set1, std::vector<int> set2, double match_weight);

	void insert(Solution sol);
	void insert(Solution sol, std::vector<int> lookup1, std::vector<int> lookup2);
	//completes the matching in the way that every non-matched vertex, a negative increasing number is assigned, s.t. they can be distinguished from successfull matches
	void completeMatching();

	int getMatchCount();

	std::pair<int, int> getMatchIDs();

	std::vector<int> getMatchIndices(bool map);

	std::vector<double> getMatchWeights(bool map);

	std::vector < std::pair<std::vector<int>, std::vector<int>>> getMatching();

	double getTargetValue();

    //returns a vector {NumUnmatched,Num1to1Matches,Num1toNMatches,NumMtoNMatches}
    std::vector<int> getMatchDistribution();

	//returns a vector {Obj1to1, Obj1toN, ObjMtoN}
	std::vector<double> getObjectiveDistribution();
};

#endif