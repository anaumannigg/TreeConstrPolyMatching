#include "../include/solution.h"

Solution::Solution() {
    this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

    //init vectors to remember which match each polygon belongs to (-1 means unmatched)
    this->set1_match_index = std::vector<int>();
    this->set1_match_weight = std::vector<double>();
    this->set2_match_index = std::vector<int>();
    this->set2_match_weight = std::vector<double>();

    this->match_count = 0;

    this->match_ids = { -1,0 };

    this->target_value = 0.0;

	this->numUnmatched = 0;
	this->num1to1 = 0;
	this->num1toN = 0;
	this->numMtoN = 0;

	this->obj1to1 = 0.0;
	this->obj1toN = 0.0;
	this->objMtoN = 0.0;

}

Solution::Solution(int num_polys1,int num_polys2) {
	this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

	//init vectors to remember which match each polygon belongs to (-1 means unmatched)
	this->set1_match_index = std::vector<int>(); this->set1_match_index.resize(num_polys1, -1);
	this->set1_match_weight = std::vector<double>(); this->set1_match_weight.resize(num_polys1, 0.0);
	this->set2_match_index = std::vector<int>(); this->set2_match_index.resize(num_polys2, -1);
	this->set2_match_weight = std::vector<double>(); this->set2_match_weight.resize(num_polys2, 0.0);

	this->match_count = 0;

	this->match_ids = { -1,0 };

	this->target_value = 0.0;

	this->numUnmatched = 0;
	this->num1to1 = 0;
	this->num1toN = 0;
	this->numMtoN = 0;

	this->obj1to1 = 0.0;
	this->obj1toN = 0.0;
	this->objMtoN = 0.0;
}

Solution::Solution(int num_polys1, int num_polys2, std::pair<int,int> match_id_limits) {
	this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

	//init vectors to remember which match each polygon belongs to (-1 means unmatched)
	this->set1_match_index = std::vector<int>(); this->set1_match_index.resize(num_polys1, -1);
	this->set1_match_weight = std::vector<double>(); this->set1_match_weight.resize(num_polys1, 0.0);
	this->set2_match_index = std::vector<int>(); this->set2_match_index.resize(num_polys2, -1);
	this->set2_match_weight = std::vector<double>(); this->set2_match_weight.resize(num_polys2, 0.0);

	this->match_count = 0;

	this->match_ids = match_id_limits;

	this->target_value = 0.0;

	this->numUnmatched = 0;
	this->num1to1 = 0;
	this->num1toN = 0;
	this->numMtoN = 0;

	this->obj1to1 = 0.0;
	this->obj1toN = 0.0;
	this->objMtoN = 0.0;
}

void Solution::addMatch(std::vector<int> set1, std::vector<int> set2, double match_weight) {
	//only add the match, if both sets are non-empty, unmatched polygons should get negative match indices via "completeMatching" in the end
	if (set1.size() == 0 || set2.size() == 0) return;

	if(match_weight < 0.0) {
		cout << "warning: adding negative match with weight " << match_weight << endl;
		cout << "SET1: ";
		for(auto& o : set1) cout << o << ", ";
		cout << "\nSET2: ";
		for(auto& a : set2) cout << a << ", ";
		cout << endl;
	}

	//remember which kind of match is added
	if (set1.size() == 1 && set2.size() == 1) {
		num1to1++; obj1to1 += match_weight;
	}
	else if (set1.size() == 1 || set2.size() == 1) {
		num1toN++; obj1toN += match_weight;
	}
	else {
		numMtoN++; objMtoN += match_weight;
	}

	std::pair<std::vector<int>, std::vector<int>> pair(set1, set2);
	this->matching.push_back(pair);
	//remember match index for each poly
	for (const auto& p : set1) {
		this->set1_match_index[p] = this->match_ids.second;
		this->set1_match_weight[p] = match_weight;
	}
	for (const auto& p : set2) {
		this->set2_match_index[p] = this->match_ids.second;
		this->set2_match_weight[p] = match_weight;
	}
	this->match_ids.second++;
	this->match_count++;

	this->target_value += match_weight;
}

void Solution::insert(Solution sol) {
	//if the solution to insert has references, we want to use the references as then we are
	//building a global solution from local solutions

	auto& matching_ids = sol.has_references ? sol.matching_references : sol.matching;

	//for each match in the argument, add the match to the solution instance
	for (int i=0; i< matching_ids.size(); i++){
		std::vector<int> set1 = matching_ids[i].first;
		std::vector<int> set2 = matching_ids[i].second;
		this->addMatch(set1, set2, sol.getMatchWeights(0)[sol.matching[i].first[0]]);
	}

}

void Solution::insert(Solution sol, const std::vector<int>& lookup1, const std::vector<int>& lookup2) {
	auto sol_matching = sol.getMatching();

	//for each match in the argument, add the match to the solution instance
	for (const auto& m : sol_matching) {
		std::vector<int> set1, set2;
		for (const auto& i : m.first) set1.push_back(lookup1[i]);
		for (const auto& i : m.second) set2.push_back(lookup2[i]);
		this->addMatch(set1, set2, sol.getMatchWeights(0)[m.first[0]]);
	}

}

void Solution::completeMatching() {
	int match_id = this->match_ids.first;
	for (int i = 0; i < this->set1_match_index.size(); i++) {
		if (this->set1_match_index[i] == -1) {
			//set match ID for unmatched polygon to negative match ID, weight does not need to be updated, as it is 0.0 as default
			this->set1_match_index[i] = match_id--;
			this->match_count++;
		}
	}
	for (int i = 0; i < this->set2_match_index.size(); i++) {
		if (this->set2_match_index[i] == -1) {
			//set match ID for unmatched polygon to negative match ID, weight does not need to be updated, as it is 0.0 as default
			this->set2_match_index[i] = match_id--;
			this->match_count++;
		}
	}
	this->numUnmatched = std::abs(match_id)-1;
	this->match_ids.first = match_id;
}

void Solution::setReferences(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2) {
	this->has_references = true;

	for (const auto& match : this->matching) {
		std::vector<int> poly_refs1,poly_refs2;
		for (const auto& p : match.first) {
			if (polys1[p].global_id == -1) std::cerr << "WARNING: global id does not exist" << std::endl;
			poly_refs1.push_back(polys1[p].global_id);
		}
		for (const auto& p : match.second) {
			if (polys2[p].global_id == -1) std::cerr << "WARNING: global id does not exist" << std::endl;
			poly_refs2.push_back(polys2[p].global_id);
		}
		this->matching_references.push_back(make_pair(poly_refs1, poly_refs2));
	}
}

int Solution::getMatchCount() {
	return this->match_count;
}

std::pair<int, int> Solution::getMatchIDs() {
	return this->match_ids;
}

//get match indices of the map indicated by the argument
std::vector<int> Solution::getMatchIndices(bool map){
	if (!map) return this->set1_match_index;
	else return this->set2_match_index;
}

std::vector<double> Solution::getMatchWeights(bool map) {
	if (!map) return this->set1_match_weight;
	else return this->set2_match_weight;
}

std::vector < std::pair<std::vector<int>, std::vector<int>>> Solution::getMatching() {
	return this->matching;
}

std::vector < std::pair<std::vector<int>, std::vector<int>>> Solution::getMatchingReferences() {
	return this->matching_references;
}

double Solution::getTargetValue() {
	return this->target_value;
}

std::vector<int> Solution::getMatchDistribution() {
    return {this->numUnmatched, this->num1to1, this->num1toN, this->numMtoN};
}

std::vector<double> Solution::getObjectiveDistribution() {
	return {this->obj1to1, this->obj1toN, this->objMtoN};
}

