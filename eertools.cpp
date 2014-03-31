/*
 * Utility library for the Entropy Evolution Rate.
 *
 *
 */
#include <unordered_map>
#include <cassert>
#include "eertools.h"

using namespace std;

bool eertools::calcComplexProb(vmat_double& complexProb, const std::string& s1, const std::string& s2){
    complexProb.clear();
	complexProb.resize(PROB_TABLE_SIZE);
	for(int i=0; i<PROB_TABLE_SIZE; i++)complexProb.at(i).resize(PROB_TABLE_SIZE, 0.0);

	int l = (s1.size() < s2.size())?(int)s1.size():(int)s2.size();

	if(s1.size() != s2.size())cerr << "calcComplexProb: wrong length." << endl;

	for(int i=0; i<l; i++)complexProb.at(s1.at(i)).at(s2.at(i)) += 1.0;

	for(int i=0; i<PROB_TABLE_SIZE; i++)for(int j=0; j<PROB_TABLE_SIZE; j++)complexProb.at(i).at(j) /= l;
	
	return true;		
}

bool eertools::calcProb(vec_double& prob, const std::string s){
	prob.clear();
	for(int i=0; i<PROB_TABLE_SIZE; i++)prob.push_back(0.0);

	for(int i=0; i<(int)s.size(); i++)prob.at(s.at(i)) += 1.0;

	for (int i=0; i<PROB_TABLE_SIZE; i++)prob.at(i) /= (double)s.size();

	return true;
}

double eertools::getEntropy(const vec_double& prob){
	double res = 0.0;

	for (int i=0; i<prob.size(); i++)
	{
		if(prob.at(i) == 0.0)continue;
		res -= prob.at(i) * log(prob.at(i));
	}
	res /= ln2;

	return max(res, 0.0);
}

double eertools::getMutualEntropy(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb){
	double res = 0.0;

	for (int i=0; i<PROB_TABLE_SIZE; i++)for (int j=0; j<PROB_TABLE_SIZE; j++)
	{
		if(prob1.at(i) == 0.0 || prob2.at(j) == 0.0 || complexProb.at(i).at(j) == 0.0)continue;
		res += complexProb.at(i).at(j) * log(complexProb.at(i).at(j)/prob1.at(i)/prob2.at(j));
	}
	res /= ln2;

	return max(res, 0.0);
}

double eertools::getEER(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb){
	double I = getMutualEntropy(prob1, prob2, complexProb);
	const double e1 = getEntropy(prob1);
	const double e2 = getEntropy(prob2);
	const double ret = 1.0 - (I/e1 + I/e2)/2.0;

	// treat nan
	if(max(max(max(e1,e2),I), 0.0) > 0.0)
	{
		return max(ret, 0.0);
	}else{
		return 0.0;
	}
}

double eertools::getEER2(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb){
	const double I = getMutualEntropy(prob1, prob2, complexProb);
	const double e1 = getEntropy(prob1);
	const double e2 = getEntropy(prob2);
	const double ret = 1.0 - I/(e1 + e2 - I);

	// treat nan
	if(max(max(max(e1,e2),I), 0.0) > 0.0)
	{
		return max(ret, 0.0);
	}else{
		return 0.0;
	}
}

double eertools::getECD(const vec_double& prob1, const vec_double& prob2, const vmat_double& complexProb){
	double I = getMutualEntropy(prob1, prob2, complexProb);

	return getEntropy(prob1) - I;
}


/*	軌道のカオス尺度を求める.
* 必要最小限の分布のみを用いた計算を行う
* 計算誤差を意識した実装
*
*	{nodes[x] | x∈[nStart,nStart+L-1]}が計算の対象
*	nodes: 軌道データ [0,1]^nodes.size
*	Pij: xn∈Ai,xn+1∈Aj xn,xn+1∈[nStart,nStart+L-1]である確率
*	Pi=Pi0+Pi1+...Pi(nStart+L-1)
*	計算コスト= O(L)+O(L)=O(L) O(MapSearch)*2*L+#{Aij}
* 必要メモリ= O(L) #{Aij}<=L
* 従来法での計算コスト= O(L)+O(M^2)=O(M^2)
* 必要メモリ= O(N^2)
*/
double eertools::getECD(const vec_double& nodes, const int& nStart, const int& L, const int& M){
	const double eps = 1.0 / M;
	unordered_map<pair<int, int>, int> Nij;		//#{x(n)∈Ai,x(n+1)∈Aj},x(n)∈[nStart,nStart+L-2]
	unordered_map<int, int> Ni;								//#{x(n)∈Ai,x(n)∈[nStart,nStart+L-2]
	unordered_map<pair<int, int>, int>::const_iterator n;
	int Ai, Aj;
	//分布Nij,Niを計算
	Aj = static_cast<int>(floor(nodes.at(nStart) / eps));
	for(int x1=1; x1<L; x1++)
	{
		Ai = Aj;
		Aj = static_cast<int>(floor(nodes.at(nStart+x1) / eps));	//j, x∈Aj Aj=[j*eps,j*eps+eps)
		Nij[pair<int, int>(Ai,Aj)]++;
		Ni[Ai]++;
	}

	//ECDを計算
	double plus = 0.0;
	double minus = 0.0;
	//double temp = 0.0;
	for(n = Nij.begin(); n != Nij.end(); n++)
	{
		const int n_ij = n->second;
		const int n_i = Ni[(n->first.first)];
		assert(n_ij != 0);
		plus += n_ij * log2(n_i);
		minus += n_ij * log2(n_ij);
		//temp += n_ij * log2(static_cast<double>(n_i)/n_ij);
	}
//const string fmt_f = "%1$-17.15f\t";
//cout << "a:" << boost::format(fmt_f) % ((plus-minus)/L) << "\t";
//cout << "n:" << boost::format(fmt_f) % (temp/L) << endl;
	return (plus-minus)/L;
}


double eertools::getECD(const string& seq){
	const int L = static_cast<int>(seq.size());
	unordered_map<pair<char, char>, int> Nij;
	unordered_map<char, int> Ni;
	unordered_map<pair<char, char>, int>::const_iterator n;
	int Ai, Aj;
	//分布Nij,Niを計算
	Aj = seq.at(0);
	for(int x1=1; x1<L; x1++)
	{
		Ai = Aj;
		Aj = seq.at(x1);
		Nij[pair<char, char>(Ai,Aj)]++;
		Ni[Ai]++;
	}

	//ECDを計算
	double plus = 0.0;
	double minus = 0.0;
	//double temp = 0.0;
	for(n = Nij.begin(); n != Nij.end(); n++)
	{
		const int n_ij = n->second;
		const int n_i = Ni[(n->first.first)];
		assert(n_ij != 0);
		plus += n_ij * log2(n_i);
		minus += n_ij * log2(n_ij);
		//temp += n_ij * log2(static_cast<double>(n_i)/n_ij);
	}
//const string fmt_f = "%1$-17.15f\t";
//cout << "a:" << boost::format(fmt_f) % ((plus-minus)/L) << "\t";
//cout << "n:" << boost::format(fmt_f) % (temp/L) << endl;
	return (plus-minus)/L;
	}


double eertools::getLYAP(const vec_double& nodes, const int& nStart, const int& L, const std::function<double (double t)> & df){
	double sum=0.0;
	for(int n=nStart; n<nStart+L; n++)sum += log(abs(df(nodes.at(n))));
	return sum / L;
}
