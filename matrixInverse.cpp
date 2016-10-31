#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <tuple>
#define WIDTH 3
const double error = 0.0000001, mean = 0.0, stddev = 1.0;
using namespace std;

class Sqmatrix
{
public:
   //Set functions
   void set( const size_t ntemp ); //Generate [n][n] matrix, initialize n
   Sqmatrix identity( );
   tuple<Sqmatrix,Sqmatrix> lu();
   Sqmatrix inverse_lu( );
   Sqmatrix inverse_gaussj( );
   bool empty( ) { if( n == 0 ) return true; else return false; } //Empty matrix condition
   //Constuctors
   Sqmatrix( ); //Default empty matrix
   explicit Sqmatrix( const size_t ntemp ) { set( ntemp ); } //Construct [n][n] matrix
   Sqmatrix( double** a, const size_t ntemp ); //Construct [n][n] matrix w/ 2d array
   Sqmatrix( double init, const size_t ntemp ); //Construct [n][n] matrix, initialize with value
   Sqmatrix( const Sqmatrix& mcopy ): Sqmatrix( mcopy.m, mcopy.n ) { }; //Copy constructor
   ~Sqmatrix(); //Destructor
   //Operators
   bool operator!=( const Sqmatrix& mat ); //Inequality operator
   Sqmatrix mul( const Sqmatrix& m2 ); //Matrix multiplication
   Sqmatrix operator*( const Sqmatrix& m2 ) { return mul( m2 ); } //Matrix multiplication
   Sqmatrix& operator=( Sqmatrix mat ); //Assignment operator
   friend ostream& operator<<( ostream& os, Sqmatrix& mat ); //Extraction operator: Print matrix

private:
   double **m;
   size_t n;

};

double** gen( size_t n )
{
   double** a = new double *[n];
   for( size_t i = 0; i < n; ++i )
      a[i] = new double [n];

   return a;
}

void del( double** a, size_t n )
{
   for( size_t i = 0; i < n; ++i )
   {  delete[] a[i];
      a[i] = nullptr;
   }
   delete[] a;
   a = nullptr;
}

double** rand_array( long long s, size_t n )
{
   default_random_engine generator( s );
   normal_distribution<double> distribution( mean, stddev );

   double** a = gen( n );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
         a[r][c] = distribution( generator );

   return a;
}

int main()
{
   unsigned int size, trials;
   cout << "Enter the size of the matrix: ";
   cin >> size;
   cout << "Enter the number of trials: ";
   cin >> trials;

   double** array;
   long long seed;
   chrono::high_resolution_clock::time_point start, end;
   chrono::microseconds timer_lu( 0 ), timer_gaussj( 0 );
   Sqmatrix inv( size ), id = inv.identity();
   unsigned int lu_dne = 0, lu_error = 0, gaussj_error = 0, singular = 0;

   for( size_t i = 1; i <= trials; ++i )
   {  if( size <= 10 )
         cout << "\nTrial " << i << "\n---------" << endl;
      else
         cout << endl << i << endl;

      seed = chrono::system_clock::now().time_since_epoch().count();

      array = rand_array( seed, size );
	  Sqmatrix mat( array, size );
      del( array, size );

      start = chrono::high_resolution_clock::now();
      inv = mat.inverse_lu();
      end = chrono::high_resolution_clock::now();
      timer_lu += chrono::duration_cast<chrono::microseconds>( end - start );

	  if( inv.empty() )
      {  ++lu_dne;
         ++lu_error;
	  }
      else if( mat * inv != id )
      {  cout << "LU inverse incorrect +1\n" << endl;
         ++lu_error;
      }
      else
         cout << "-\n" << endl; //Passed

      start = chrono::high_resolution_clock::now();
      inv = mat.inverse_gaussj();
      end = chrono::high_resolution_clock::now();
      timer_gaussj += chrono::duration_cast<chrono::microseconds>( end - start );

      if( inv.empty() )
      {  ++singular;
         ++gaussj_error;
	  }
	  else if( mat * inv != id )
      {  cout << "Gauss Jordan incorrect +1" << endl;
         ++gaussj_error;
      }
      else
	     cout << "-" << endl; //Passed
   }

   cout << "\nInverse LU: ";
   if( timer_lu.count() < 1000 )
      cout << timer_lu.count() << " microseconds" << endl;
   else if( timer_lu.count() >= 1000 and timer_lu.count() < 1000000 )
      cout << chrono::duration_cast<chrono::milliseconds>(timer_lu).count() << " milliseconds" << endl;
   else if( timer_lu.count() >= 1000000 and timer_lu.count() < 60000000 )
      cout << chrono::duration_cast<chrono::seconds>(timer_lu).count() << " seconds" << endl;
   else
      cout << chrono::duration_cast<chrono::minutes>(timer_lu).count() << " minutes" << endl;

   cout << "LU dne " << lu_dne << "\ninv error " << lu_error << endl;

   cout << "\nGauss Jordan: ";
   if( timer_gaussj.count() < 1000 )
      cout << timer_gaussj.count() << " microseconds" << endl;
   else if( timer_gaussj.count() >= 1000 and timer_gaussj.count() < 1000000 )
      cout << chrono::duration_cast<chrono::milliseconds>(timer_gaussj).count() << " milliseconds" << endl;
   else if( timer_gaussj.count() >= 1000000 and timer_gaussj.count() < 60000000 )
      cout << chrono::duration_cast<chrono::seconds>(timer_gaussj).count() << " seconds" << endl;
   else
      cout << chrono::duration_cast<chrono::minutes>(timer_lu).count() << " minutes" << endl;

   cout << "singular " << singular << "\ninv error " << gaussj_error << endl;

   return 0;
}

void Sqmatrix::set( const size_t ntemp )
{
   n = ntemp;
   m = new double *[n];
   for( size_t i = 0; i < n; ++i )
      m[i] = new double [n];

}

Sqmatrix Sqmatrix::identity( )
{
   Sqmatrix id( n );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
      {  if( c == r )
            id.m[r][c] = 1;
         else
            id.m[r][c] = 0;
      }
   
   return id;
}

tuple<Sqmatrix,Sqmatrix> Sqmatrix::lu( ) //LU decomposition using Gaussian elimination: permutation matrix to be implemented
{
   Sqmatrix lt( n );
   Sqmatrix ut = *this;
   double factor;

   for( size_t i = 0; i < n; ++i )
   {  if( abs(ut.m[i][i]) < error ) //Nonexistent LU condition
      {  cout << "LU dne\n" << endl;
         return forward_as_tuple( Sqmatrix(), Sqmatrix() );
      }
      lt.m[i][i] = 1;
      for( size_t j = i + 1; j < n; ++j )
      {  if( j < n )
            lt.m[i][j] = 0;
         factor = ut.m[j][i] / ut.m[i][i];
         lt.m[j][i] = factor;
         if( abs(ut.m[j][i]) > error )
            for( size_t k = 0; k < n; ++k )
                ut.m[j][k] -= ut.m[i][k] * factor;
         
      }
   }
   return forward_as_tuple( lt, ut ); //C++11 tuple return
}

Sqmatrix Sqmatrix::inverse_lu( )
{
   Sqmatrix inv = identity();
   Sqmatrix lt( n ), ut( n );
   tie( lt, ut ) = lu();

   if( lt.n == 0 and ut.n == 0 )
      return Sqmatrix();

   for( size_t i = 0; i < n; ++i )
   {  for( size_t j = 0; j < n; ++j ) //Forwards substitution: L * (L^-1) = I
      {  for( size_t k = 0; k < j + 1; ++k )
            if( k != j )
               inv.m[j][i] -= lt.m[j][k] * inv.m[k][i];
         inv.m[j][i] /= lt.m[j][j];
      }
      for( int j = n - 1; j > -1; --j ) //Backwards substitution: U * (A^-1) = L^-1
      {  for( int k = n - 1; k > j - 1; --k )
            if( k != j )
               inv.m[j][i] -= ut.m[j][k] * inv.m[k][i];
         inv.m[j][i] /= ut.m[j][j];
      }
   }
   return inv;
}

Sqmatrix Sqmatrix::inverse_gaussj( ) //Augmented matrix (implemented w/ full scaled pivoting)
{
   Sqmatrix inv = identity();
   Sqmatrix mat = *this;
   size_t max_r, max_c;
   double factor, scale[n-1], temp;

   for( size_t i = 0; i < n; ++i )
   {  /*max_c = i;
      for( size_t j = i; j < n; ++j ) //Scale comparisons
      {  max_r = i;
         for( size_t k = i + 1; k < n; ++k )
            if( abs(mat.m[j][max_r]) < abs(mat.m[j][k]) )
               max_r = k;
         scale[j] = mat.m[j][i] / mat.m[j][max_r];

         if( max_c != j and abs(scale[max_c]) < abs(scale[j]) )
            max_c = j;
      }
      
      if( abs(mat.m[max_c][i]) < error ) //Singular matrix condition
      {  cout << "Singular matrix\n" << endl;
         return Sqmatrix();
      }

      if( i != max_c ) //Row swaps
         for( size_t j = 0; j < n; ++j )
         {  temp = mat.m[i][j];
            mat.m[i][j] = mat.m[max_c][j];
            mat.m[max_c][j] = temp;
         
            temp = inv.m[i][j];
            inv.m[i][j] = inv.m[max_c][j];
            inv.m[max_c][j] = temp;
         }*/


	  if( abs(mat.m[i][i]) < error ) //Pivoting
         for( size_t j = i + 1; j < n; ++j )
         {  if( abs(mat.m[j][i]) < error )
               continue;
            factor = ( mat.m[i][i] - 1 ) / mat.m[j][i];
            for( size_t k = i; k < n; ++k )
               mat.m[i][k] -= mat.m[j][k] * factor;
            for( size_t k = 0; k < n; ++k )
			   inv.m[i][k] -= inv.m[j][k] * factor;
            break;
         }

      if( abs(mat.m[i][i]) < error ) //Singular matrix condition
      {  cout << "Singular matrix\n" << endl;
         return Sqmatrix();
      }


      if( mat.m[i][i] < 1.0 - error or mat.m[i][i] > 1.0 + error ) //Pivoting
      {  factor = mat.m[i][i];
         for( size_t j = 0; j < n; ++j )
         {  mat.m[i][j] /= factor;
            inv.m[i][j] /= factor;
         }
      }
      for( size_t j = 0; j < n; ++j ) //Eliminating
      {  if( j == i or abs(mat.m[j][i]) < error )
            continue;
         factor = mat.m[j][i];
         for( size_t k = i; k < n; ++k )
            mat.m[j][k] -= mat.m[i][k] * factor;
         for( size_t k = 0; k < n; ++k )
            inv.m[j][k] -= inv.m[i][k] * factor;
      }
   }
   return inv;
}

Sqmatrix::Sqmatrix( ) //Default declaration
{
   set( 0 );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
         m[r][c] = 0;

}

Sqmatrix::Sqmatrix( double** a, const size_t ntemp )
{
   set( ntemp );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
         m[r][c] = a[r][c]; 

}

Sqmatrix::Sqmatrix( double init, const size_t ntemp )
{
   set( ntemp );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
         m[r][c] = init;

}

Sqmatrix::~Sqmatrix( )
{
   for( size_t i = 0; i < n; ++i )
   {  delete[] m[i];
      m[i] = nullptr;
   }
   delete[] m;
   m = nullptr;
   n = 0;

}

bool Sqmatrix::operator!=( const Sqmatrix& mat )
{
   if( n != mat.n )
      return true;
   for( size_t i = 0; i < n; ++i )
      for( size_t j = 0; j < n; ++j )
         if( abs(m[i][j] - mat.m[i][j]) > error )
            return true;

   return false;
}

Sqmatrix Sqmatrix::mul( const Sqmatrix& m2 )
{
   Sqmatrix prod( n );
   for( size_t r = 0; r < n; ++r )
      for( size_t c = 0; c < n; ++c )
      {  prod.m[r][c] = 0;
         for( size_t k = 0; k < n; ++k )
            prod.m[r][c] += m[r][k] * m2.m[k][c];
      }

   return prod;
}

Sqmatrix& Sqmatrix::operator=( Sqmatrix mat )
{  //Swap n
   size_t ntemp = n;
   n = mat.n;
   mat.n = ntemp;

   //Swap m
   double** mtemp = m;
   m = mat.m;
   mat.m = mtemp;

   return *this;
}

ostream& operator<<( ostream& os, Sqmatrix& mat )
{
   for( size_t r = 0; r < mat.n; ++r )
   {  os << "[ ";
      for( size_t c = 0; c < mat.n; ++c )
      {  if( abs(mat.m[r][c]) < error )
            os << setw(WIDTH) << "0" << "  ";
         else
            os << setw(WIDTH) << setprecision(3) << mat.m[r][c] << "  ";
      }
      os << " ]" << endl;
   }
   return os;
}
