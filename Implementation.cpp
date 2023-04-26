// TOOD: Redo with a Q object so that it stores the lambda internally 

#include <vector>
#include <cmath>
#include <random>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <map>

#include <Rcpp.h>

constexpr double infinity = std::numeric_limits<double>::infinity();

const int N = 100;

const double MIN_IB = 0.0;
const double MAX_IB = 10.0;

// for Golden-section search
const auto phi = (1+sqrt(5.))/2.;
const auto phi_inv = 1.0/phi;

// [[Rcpp::export]] 
double logplusexp(double a, double b) {
	double mx = std::max(a,b);
// 	return mx + log(exp(a-mx)+exp(b-mx));
	return mx + log1p(exp(std::min(a,b)-mx));
}

// [[Rcpp::export]] 
double logminusexp(double a, double b) {
	// log(exp(a)-exp(b)) 
	// log(1-exp(b)/exp(a))
	return a + log1p(-exp(b-a));
}

// [[Rcpp::export]] 
double logsumexp(Rcpp::NumericVector& v) {
	double m = Rcpp::max(v);
	double s = 0.0;
	for(auto& x : v) { s += exp(x-m); }
	return m + log(s);
}

double rnd(double x, int n) {
	int m = pow(10,n);
	return round(x*m)/m;
}

// our prior over numbers is proportional to 1/n^2
// NOTE: this is not normalized over 1..N, but up to infinity
// [[Rcpp::export]] 
double P(int n) {
	assert(n >= 1);
	return pow(n,-2) / 1.6449;
}

// for a given lambda_n (lambda depends on n) this returns
// the *log* probability of responding k to n, *unnormalized*
// This is (5) in Cheyette & Piantadosi
// [[Rcpp::export]] 
double logQ(int k, int n, double lambda_n) {
	assert(k >= 1);
	assert(n >= 1);
	assert(lambda_n > 0);
	return log(P(k)) - P(n) * pow(n-k,2) / lambda_n;
}

// what is the log normalizer for the Q distribution
// again treating N as the largest possible 
// [[Rcpp::export]] 
double logQz(int n, double lambda_n) {
	assert(n >= 1);
	assert(lambda_n > 0);
	double z = -infinity;
	for(int k=1;k<=N;k++) {
		z = logplusexp(z, logQ(k,n,lambda_n));
	}
	return z;
}

// find the KL divergence between Q and the prior
// summing up to N
// [[Rcpp::export]] 
double KL(int n, double lambda_n) {
	assert(n >= 1);
	assert(lambda_n > 0);
	double qz = logQz(n, lambda_n); // cache this
	double kl = 0.0;
	for(int k=1;k<=N;k++) {
		double lq = logQ(k,n,lambda_n) - qz; // Q(k|n)
		kl += exp(lq) * (lq - log(P(k)));
	}
	return kl;
}



// convert an information bound ib into a lambda_n 
// considering the numbers going up to N
// for this we use a binary search, which works because our
// function of lambda should be decreasing
// [[Rcpp::export]] 
double find_lambda_n(const double ib, const int n) {
	assert(n >= 1); assert(ib >= 0.0);
	const double tolerance=0.000001;
	double lower=tolerance; // if this is too large, everything breaks
	double upper=100000.0; // if this is too small, everything breaks
	double midpoint{0.0}; 
	/*	
	for(double l=0.0001;l<10000;l *= 2.) {
		std::cout << l <<"\t"<< KL(n,l) << std::endl;
	}
	*/
	// first check, we may not need to search if
	// the information bound is high enough that
	// the KL cant be that big
// 	if(KL(n,upper) > ib) return upper;
// 	if(KL(n,lower) < ib) return lower;
	
	// FOR NOW: Fixed number of iterations
	// we can expect a precision like (uppper-lower) / 2^iterations
	while(upper-lower > tolerance){
		assert(upper > lower);
		midpoint = lower+(upper-lower)/2;
		
		auto fmid = KL(n,midpoint);
		//std::cout << "FMID=" << fmid <<"\t"<< lower <<"\t"<< upper <<"\t"<< ib << std::endl;
		
		// the KL divergence should be (generally) a decreasing function of lambda
		if(fmid > ib) { lower = midpoint; }
		else		  { upper = midpoint; }
	}
	
	return midpoint;
}

// Represent a Q distribution centered on n
// NOTE that initializing this is not trivial because it involves finding a lambda and a normalizer 
struct QDistribution {
	
	const int n;
	double lambda;
	double z;
	const double ib;
	
	QDistribution(int _n, double _ib) : n(_n), ib(_ib) {
		lambda = find_lambda_n(ib,n);
		z = logQz(n,lambda);
	}
	
	double lp(const int k) const {
		return logQ(k,n,lambda)-z;
	}
	
	void show(int m) const {
		for(int k=1;k<m;k++) {
			std::cout << n <<":\t"<< k <<"\t"<< exp(lp(k)) <<"\t"<< lambda << std::endl;
		}
	}
	
	int sample() const {
		auto rnd = R::runif(0,1);
		for(size_t k=1;k<N;k++) {
			auto p = exp(lp(k));
			rnd = rnd - p; 
			if(rnd < 0)
				return k;
		}
		assert(0);
		return 0;
	}
	
};

// [[Rcpp::export]] 
double discrimination_likelihood(const Rcpp::NumericVector& n1, 
								 const Rcpp::NumericVector& n2, 
								 const Rcpp::NumericVector& accuracy, 
								 const double ib) {
	assert(n1.size() == n2.size());
	assert(n1.size() == accuracy.size());
	
	// TODO: Cache lambda_n
	double ll = 0.0;
	for(size_t i=0;i<n1.size();i++) {
		
		const QDistribution q1(n1.at(i),ib);
		const QDistribution q2(n2.at(i),ib);
		
		// add up all the probability that you think n1 > n2
		// which is a sum over what n1 is, times the prob that n2 is less
		auto lpn1gtn2 = -infinity; 
		
		// first add up the probability that we represent the same numbers
		// and that's assumed to be 50% accuracy 
		for(int r=1;r<N;r++) {
			lpn1gtn2 = logplusexp(lpn1gtn2, q1.lp(r)+q2.lp(r)+log(0.5));
		}
		
		// now add up the probability of getting any r2 < r1, starting with r2=1
		auto r2cdf = q2.lp(1);
		for(int r=2;r<N;r++) {
			//std::cout << ib <<"\t"<< n1.at(i) << "\t" << n2.at(i) << "\t" << r  <<"\t"<< r2cdf <<"\t"<< lpr1 <<"\t"<< 	
			//			 lpn1gtn2 <<"\t"<< accuracy.at(i) << std::endl;
			
			// add up P(r|n1,n2) P(r2<r|n1,n2)
			lpn1gtn2 = logplusexp(lpn1gtn2, q1.lp(r) + r2cdf );
			
			// update P(r2 < r) by adding up the probability
			r2cdf = logplusexp(r2cdf, q2.lp(r));
		}
		
		if(accuracy.at(i) == (n1.at(i) > n2.at(i))) {
			ll += lpn1gtn2;
		}
		else {
			ll += log1p(-exp(lpn1gtn2)); //log(1-exp(lp));
		}
		//std::cout << ib <<"\t"<< lpn1gtn2 <<"\t"<< accuracy.at(i) <<"\t"<< ll << std::endl;
	}
// 	std::cout << ib <<"\t"<< ll <<"\t"<< n1.size() << std::endl;
	return ll;
}

// [[Rcpp::export]] 
Rcpp::NumericVector estimation_sample(const Rcpp::NumericVector& shown, 
									  const double ib) {
	assert(ib >= 0.0);
	
	Rcpp::NumericVector out;
	for(size_t i=0;i<shown.size();i++) {
		const QDistribution q(shown.at(i), ib);
	
		// multinomial sample 
		auto rnd = R::runif(0,1);
		for(size_t k=1;k<N;k++) {
			auto p = exp(q.lp(k));
// 			std::cout << k <<"\t"<< p <<"\t" << rnd <<"\t" << std::endl;
			rnd = rnd - p; 
			if(rnd < 0) {
				out.push_back(k);
				break;
			}
		}
		assert(0); // should not get here
	}
	return out;
}

// [[Rcpp::export]] 
Rcpp::NumericVector discrimination_sample(const Rcpp::NumericVector& n1,
									      const Rcpp::NumericVector& n2, 
										  const double ib) {
	// returns estimates of *accuracy* in discrimination 
	assert(n1.size() == n2.size());
	assert(ib >= 0.0);
	
	// conveniently we can just use likelihood for this 
	Rcpp::NumericVector out; 
	for(size_t i=0;i<n1.size();i++) {
		assert(n1.at(i) != n2.at(i)); // 
		const QDistribution q1(n1.at(i),ib);
		const QDistribution q2(n2.at(i),ib);
		out.push_back( (q1.sample() > q2.sample()) == (n1.at(i) > n2.at(i)));
	}
	return out;
}



// [[Rcpp::export]] 
double estimation_likelihood(const Rcpp::NumericVector& shown, 
							 const Rcpp::NumericVector& response,
							 const double ib) {
	assert(response.size() == shown.size());
	assert(ib >= 0.0);
	
	// TODO: Cache lambda_n 
	double ll = 0.0;
	for(size_t i=0;i<response.size();i++) {
		const QDistribution q(shown.at(i), ib);
		ll += q.lp(response[i]);
	}
	return ll;
}

// [[Rcpp::export]]
double estimation_fit_ib_bound(const Rcpp::NumericVector& shown,
							   const Rcpp::NumericVector& response) {
	// Implement a golden-section search here to fit. 
	// this searches in [a,b] with interior points a < x < y < b
	// https://en.wikipedia.org/wiki/Golden-section_search

	const auto tolerance=0.001;
	
	auto a  = MIN_IB;  
	auto b  = MAX_IB; 
	
	while(b-a > tolerance){
		assert(b > a);
// 		std::cout << a << "\t" << b << std::endl;
		
		const auto x = b - phi_inv*(b-a);
		const auto y = a + phi_inv*(b-a);
		
		auto fx = estimation_likelihood(shown,response,x);
		auto fy = estimation_likelihood(shown,response,y);
		
		if( std::max(fx,fy)/std::min(fx,fy) < 1.01 and fy < response.size()*(-10) ) { 
			// first a special case here -- if fx \approx fy and the estimation_likelihood is low,
			// then we probably have an ib that is too high, so treat this as b is bad
			b = y;
		}
		else if(fx > fy) { 
			b = y; 
		}
		else { 
			a = x; 
		}
	}
	
	return a + (b-a)/2.;
}



// [[Rcpp::export]]
double discrimination_fit_ib_bound(const Rcpp::NumericVector& n1, 
								   const Rcpp::NumericVector& n2,
								   const Rcpp::NumericVector& accuracy) {
	// Implement a golden-section search here to fit. 
	// this searches in [a,b] with interior points a < x < y < b
	// https://en.wikipedia.org/wiki/Golden-section_search

	const auto tolerance=0.001;
	
	auto a  = MIN_IB;  
	auto b  = MAX_IB; 
	
	while(b-a > tolerance){
		assert(b > a);
		
		const auto x = b - phi_inv*(b-a);
		const auto y = a + phi_inv*(b-a);
		assert(x < y);
		auto fx = discrimination_likelihood(n1,n2,accuracy,x);
		auto fy = discrimination_likelihood(n1,n2,accuracy,y);
		//std::cout << a << "\t" << b <<"\t"<< fx <<"\t"<< fy << std::endl;
		
		if(fy <= -infinity or std::isnan(fy)) {
			b = y;
		}
		else if(fx <= -infinity or std::isnan(fx)) {
			a = x;
		}
		else if(fx > fy) { 
			b = y; 
		}
		else { 
			a = x; 
		}
	}
	
	return a + (b-a)/2.;
}



