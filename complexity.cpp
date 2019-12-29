/* Copyright 2019 K. J. Teng

 Permission is hereby granted, free of charge,
 to any person obtaining a copy of this software
 and associated documentation files (the "Software"),
 to deal in the Software without restriction,
 including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit
 persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice
 shall be included in all copies or substantial portions
 of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
 ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
 TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
 SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
 OR OTHER DEALINGS IN THE SOFTWARE.

 */

/* Iterator wrapper */

#include <iostream>
#include <iterator>

struct IteratorData {
	size_t access_count, move_count, compare_count;
	IteratorData(): access_count(0), move_count(0), compare_count(0) {}
	void reset() { access_count = move_count = compare_count = 0; }
};

template <class T>
class MyIterator {
public:
	typedef typename std::iterator_traits<T>::difference_type   difference_type;
	typedef typename std::iterator_traits<T>::value_type        value_type;
	typedef typename std::iterator_traits<T>::pointer           pointer;
	typedef typename std::iterator_traits<T>::reference         reference;
	typedef typename std::iterator_traits<T>::iterator_category iterator_category;
	
	MyIterator(T it, IteratorData *dt): iterator(it), data(dt) {}
	
	/* Input iterator / Forward iterator */
	reference   operator * () { data->access_count++; return *iterator; }
	pointer     operator ->() { data->access_count++; return iterator.operator->(); }
	bool        operator ==(const MyIterator &other) const { data->compare_count++; return iterator == other.iterator; }
	bool        operator !=(const MyIterator &other) const { data->compare_count++; return iterator != other.iterator; }
	MyIterator &operator ++() { data->move_count++; iterator++; return *this; }
	MyIterator  operator ++(int) { data->move_count++; return MyIterator(iterator++, data); }
	
	/* Bidirectional iterator */
	MyIterator &operator --() { data->move_count++; iterator--; return *this; }
	MyIterator  operator --(int) { data->move_count++; return MyIterator(iterator--, data); }
	
	/* Random-access iterator */
	MyIterator &operator +=(difference_type n) { data->move_count++; iterator += n; return *this; }
	MyIterator  operator + (difference_type n) const { data->move_count++; return MyIterator(iterator + n, data); }
	MyIterator &operator -=(difference_type n) { data->move_count++; iterator -= n; return *this; }
	MyIterator  operator - (difference_type n) const { data->move_count++; return MyIterator(iterator - n, data); }
	difference_type operator - (const MyIterator &other) const { data->compare_count++; return iterator - other.iterator; }
	reference   operator [](difference_type n) { data->move_count++; data->access_count++; return iterator[n]; }
	bool        operator < (const MyIterator &other) const { data->compare_count++; return iterator < other.iterator; }
	bool        operator > (const MyIterator &other) const { data->compare_count++; return iterator > other.iterator; }
	bool        operator <=(const MyIterator &other) const { data->compare_count++; return iterator <=other.iterator; }
	bool        operator >=(const MyIterator &other) const { data->compare_count++; return iterator >=other.iterator; }
private:
	T iterator;
	IteratorData *data;
};

template <class T>
MyIterator<T> operator + (typename MyIterator<T>::difference_type n, MyIterator<T> it) { return it + n; }

template <class T>
MyIterator<T> wrapIterator(T orig, IteratorData &data) {
	return MyIterator<T>(orig, &data);
}

/* Statistics, using least-squares method */

#include <cmath>

template <class T>
class Statistics {
public:
	Statistics(): m_weight(0), m_sum_x(0), m_sum_sq_x(0), m_sum_y(0), m_sum_sq_y(0), m_sum_xy(0) {}
	
	void addData(T x, T y, T weight = 1) {
		T weightedX = x * weight, weightedY = y * weight;
		m_weight += weight;
		m_sum_x += weightedX;
		m_sum_sq_x += weightedX * x;
		m_sum_y += weightedY;
		m_sum_sq_y += weightedY * y;
		m_sum_xy += weightedX * y;
	}
	
	void removeData(T x, T y, T weight = 1) {
		addData(x, y, -weight);
	}
	
	void clear() {
		m_weight = m_sum_x = m_sum_sq_x = m_sum_y = m_sum_sq_y = m_sum_xy = 0;
	}
	
	T getN() const { return m_weight; }
	T getSumX() const { return m_sum_x; }
	T getSumSqX() const { return m_sum_sq_x; }
	T getSumY() const { return m_sum_y; }
	T getSumSqY() const { return m_sum_sq_y; }
	T getSumXY() const { return m_sum_xy; }
	
	T getAvgX() const { return m_sum_x / m_weight; }
	T getAvgY() const { return m_sum_y / m_weight; }
	T getVarX() const { return m_sum_sq_x - m_sum_x * m_sum_x / m_weight; }
	T getVarY() const { return m_sum_sq_y - m_sum_y * m_sum_y / m_weight; }
	T getVarXY() const { return m_sum_xy - m_sum_x * m_sum_y / m_weight; }
	T getDevX() const { return std::sqrt(getVarX()); }
	T getDevY() const { return std::sqrt(getVarY()); }
	T getDevXY() const { return std::sqrt(getVarXY()); }
	
	T getR() const {
		return (m_weight * m_sum_xy - m_sum_x * m_sum_y) /
		       std::sqrt((m_weight * m_sum_sq_x - m_sum_x * m_sum_x) * (m_weight * m_sum_sq_y - m_sum_y * m_sum_y));
	}
	T getB() const { return (m_weight * m_sum_xy - m_sum_x * m_sum_y) / (m_weight * m_sum_sq_x - m_sum_x * m_sum_x); }
	T getA() const { return (m_sum_y - getB() * m_sum_x) / m_weight; }
	
private:
	T m_weight, m_sum_x, m_sum_sq_x, m_sum_y, m_sum_sq_y, m_sum_xy;
};

/* Complexity models. You can add other models below. */

#include <iostream>

class ComplexityModel {
public:
	friend std::ostream & operator << (std::ostream &, const ComplexityModel &);
	virtual void addData(double n, double Fn) = 0;
	virtual ~ComplexityModel() {}
	double getR() const { return m_stats.getR(); }
	
	bool operator < (const ComplexityModel &other) const { return getR() < other.getR(); }
protected:
	ComplexityModel() {}
	virtual void print(std::ostream &) const = 0;
	Statistics<double> m_stats;
};

std::ostream & operator << (std::ostream &os, const ComplexityModel &model) { model.print(os); return os; }

class LogarithmComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(std::log(n), Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " ln(n)"; }
};

class SqrtComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(std::sqrt(n), Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " sqrt(n)"; }
};

class LinearComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n, Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n"; }
};

class NLogNComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n * std::log(n), Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n ln(n)"; }
};

class NLog2NComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n * std::log(n) * std::log(n), Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n ln^2(n)"; }
};

class N1_5Complexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(std::pow(n, 1.5), Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n^(3/2)"; }
};

class SquareComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n * n, Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n^2"; }
};

class CubeComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n * n * n, Fn); }
	void print(std::ostream &os) const { os << m_stats.getB() << " n^3"; }
};

class ExponentialComplexity: public ComplexityModel {
	void addData(double n, double Fn) { m_stats.addData(n, std::log(Fn)); }
	void print(std::ostream &os) const { os << std::exp(m_stats.getA()) << '*' << std::exp(m_stats.getB()) << "^n"; }
};

template <class T>
struct ptr_compare {
	bool operator () (const T * lhs, const T * rhs) const { return *lhs < *rhs; }
};

#include <algorithm>

#define MODEL_NUM 9
class ComplexityAnalysis {
	ComplexityModel *models[4][MODEL_NUM];
public:
	ComplexityAnalysis() {
		for (int i = 0; i < 4; i++) {
			models[i][0] = new LogarithmComplexity;
			models[i][1] = new SqrtComplexity;
			models[i][2] = new LinearComplexity;
			models[i][3] = new NLogNComplexity;
			models[i][4] = new NLog2NComplexity;
			models[i][5] = new N1_5Complexity;
			models[i][6] = new SquareComplexity;
			models[i][7] = new CubeComplexity;
			models[i][8] = new ExponentialComplexity;
		}
	}
	
	void addData(size_t n, const IteratorData &data) {
		for (int i = 0; i < MODEL_NUM; i++) {
			models[0][i]->addData(n, data.access_count);
			models[1][i]->addData(n, data.move_count);
			models[2][i]->addData(n, data.compare_count);
			models[3][i]->addData(n, data.access_count + data.move_count + data.compare_count);
		}
	}
	
	void getResult() {
		ComplexityModel *bestFit[4] = {nullptr};
		for (int i = 0; i < 4; i++) {
			bestFit[i] = *std::max_element(models[i] + 0, models[i] + MODEL_NUM, ptr_compare<ComplexityModel>());
		}
		std::cout << "Access count best fit: " << *bestFit[0] << ", R = " << bestFit[0]->getR() << std::endl;
		std::cout << "Move count best fit: " << *bestFit[1] << ", R = " << bestFit[1]->getR() << std::endl;
		std::cout << "Compare count best fit: " << *bestFit[2] << ", R = " << bestFit[2]->getR() << std::endl;
		std::cout << "Overall best fit: " << *bestFit[3] << ", R = " << bestFit[3]->getR() << std::endl;
	}
	
	~ComplexityAnalysis() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < MODEL_NUM; j++) {
				delete models[i][j];
			}
		}
	}
};


/* Test code -- test std::sort's complexity */

#include <algorithm>
#include <list>
#include <vector>
#include <random>
#include <ctime>

#define MAX ((1 << 22) + 1)

using namespace std;

struct MySequence {
	MySequence(int x = 0): num(x), state(-1) {}
	int operator() () { return num += state; }
	int num, state;
};

int main() {
	size_t seed = time(NULL);
	cerr << "seed = " << seed << endl;
	minstd_rand generator;
	
	vector<int> assist(MAX);
	
	int seq = 1;
	ComplexityAnalysis ca;
	for (int sz = 4; sz < MAX; sz = (sz << 1)) {
		std::vector<int> container;
		cout << "** Case " << seq++ << ": " << endl;
		cout << "size = " << sz << endl;
		// randomly generate a sequence with sz numbers
		container.resize(sz);
		generate(container.begin(), container.end(), generator);
		
		IteratorData data; // test data
		// begin and end are iterators
		auto begin = wrapIterator(container.begin(), data), end = wrapIterator(container.end(), data);
		// record start time
		clock_t start_time = clock();
		sort(begin, end, wrapIterator(assist.begin(), data));
		// record end time
		clock_t end_time = clock();
		
		// check if the sequence is sorted
		auto until = is_sorted_until(container.begin(), container.end());
		if (until != container.end()) {
			cout << "Sorting failed -- sorted len = " << until - container.begin() << endl;
			continue;
		}
		
		// output single-case result
		cout << "# access = " << data.access_count << endl;
		cout << "# compare = " << data.compare_count << endl;
		cout << "# move = " << data.move_count << endl;
		cout << "time = " << (end_time - start_time) / (double)CLOCKS_PER_SEC << " second(s)" << endl;
		cout << endl;
		
		// add to ComplexityAnalysis to fit a model
		ca.addData(sz, data);
	}
	
	// output the best-fit model
	ca.getResult();
	return 0;
}
