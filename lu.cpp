#include<iostream>
#include<vector>
#include<iomanip>
using namespace std;

void cout_matrix(vector<vector<double>> A){
	int m=A.size();  //редове
	int n=A.at(0).size(); //стълбове
	for (int i=0;i<m;i++){
		for (int j=0;j<n;j++){
			cout<<setprecision(14)<<A.at(i).at(j)<<" ";
		}
		cout<<endl;
	}
}

void cout_vector(vector<double> v){
	int n=v.size();
	for (int i=0;i<n;i++){
		cout<<v.at(i)<<endl;
	}
}

vector<vector<double>> matmul(vector<vector<double>> A,vector<vector<double>> B){
	int m=A.size();
	int k=B.size();
	int n=B.at(0).size();
	vector<vector<double>> C; //C=AB, c_ij=a_it b_tj
	for (int i=0;i<m;i++){
		vector<double> rowi;
		for (int j=0;j<n;j++){
			double cij=0;
			for (int t=0;t<k;t++){
				cij+=A.at(i).at(t)*B.at(t).at(j);
			}
			rowi.push_back(cij);
		}
		C.push_back(rowi);
	}
	return C;
}

vector<vector<double>> zero_matrix(int n){
	vector<vector<double>> res;
	vector<double> row;
	for (int j=0;j<n;j++){
		row.push_back(0);
	}
	for (int i=0;i<n;i++){
		res.push_back(row);
	}
	return res;
}

vector<double> zero_vector(int n){
	vector<double> res;
	for (int i=0;i<n;i++){
		res.push_back(0);
	}
	return res;
}

vector<vector<double>> identity(int n){
	vector<vector<double>> res;
	for (int i=0;i<n;i++){
		vector<double> rowi;
		for (int j=0;j<n;j++){
			if (i==j){rowi.push_back(1);}
			else {rowi.push_back(0);}
		}
		res.push_back(rowi);
	}
	return res;
}

double  _sum(int i,int j, int p, vector<vector<double>> L,vector<vector<double>> U){
	double res=0;
	for (int k=0;k<=p;k++){
		res+=L.at(i).at(k)*U.at(k).at(j);
	}
	return res;
}

void lu(vector<vector<double>> A,vector<vector<double>> &L,vector<vector<double>> &U){
	int n=A.size();
	L=zero_matrix(n);
	U=identity(n);
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (j==0){
				L.at(i).at(j)=A.at(i).at(j);
			}
			else if (i>=j){
				L.at(i).at(j)=A.at(i).at(j)-_sum(i,j,j-1,L,U);
			}
			else {
				U.at(i).at(j)=(1/L.at(i).at(i))*(A.at(i).at(j)-_sum(i,j,i-1,L,U));
			}
		}
	}
}

double det(vector<vector<double>> A){
	int n=A.size();
	vector<vector<double>> L,U;
	lu(A,L,U);
	double res=1;
	for (int k=0;k<n;k++){
		res*=L.at(k).at(k);
	}
	return res;
}

vector<double> solve(vector<vector<double>> A,vector<double> b){
	int n=b.size();
	vector<double> x=zero_vector(n);
	vector<double> y=zero_vector(n);
	vector<vector<double>> L,U;
	lu(A,L,U);
	for (int j=0;j<n;j++){
		double yj=b.at(j);
		for (int k=0;k<=j-1;k++){
			yj-=L.at(j).at(k)*y.at(k);
		}
		y.at(j)=yj/L.at(j).at(j);
	}
	for (int i=n-1;i>=0;i--){
		double xi=y.at(i);
		for (int k=i+1;k<n;k++){
			xi-=U.at(i).at(k)*x.at(k);
		}
		x.at(i)=xi;
	}
	return x;
}

int main(){
	vector<vector<double>> L,U;
	vector<vector<double>> A={{6,18,3},{2,12,1},{4,15,3}};
	vector<double> b={1,2,3};
	lu(A,L,U);
	cout<<"L:"<<endl;
	cout_matrix(L);
	cout<<endl;
	cout<<"U:"<<endl;
	cout_matrix(U);
	cout<<"det A = "<<det(A)<<endl;
	cout<<"LU:"<<setprecision(10)<<endl;
	cout_matrix(matmul(L,U));
	cout<<"Ax=b, x=:"<<endl;
	vector<double> x=solve(A,b);
	cout_vector(solve(A,b));
	return 0;
}
