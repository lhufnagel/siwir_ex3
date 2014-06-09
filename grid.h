namespace lbm
{

	template< typename Type, uint Cellsize >
		class Grid
		{
			public:
				inline Grid();
				inline Grid( uint xsize, uint ysize );
				inline ~Grid();

				inline Type& operator()( uint x, uint y, uint f );
				inline Type  operator()( uint x, uint y, uint f ) const;

				void swap(Grid<Type,Cellsize>& grid) /* throw() */
				{
					std::swap( xsize_, grid.xsize_ );
					std::swap( ysize_, grid.ysize_ );
					std::swap( data_, grid.data_ );
					}
			private:
				uint xsize_;  // Number of nodes in x-dimension
				uint ysize_;  // Number of nodes in y-dimension
				Type* data_;  // Linearized, 1-dimensional representation of the 2D data grid
		};


	// Implementation of the default constructor
	template< typename Type, uint Cellsize > 
		Grid<Type,Cellsize>::Grid():
			xsize_(0),
			ysize_(0),
			data_(0)
	{}

	// Implementation of the initialization constructor
	template< typename Type, uint Cellsize >
		Grid<Type,Cellsize>::Grid( uint xsize, uint ysize ):
			xsize_(xsize),
			ysize_(ysize),
			data_(new Type[Cellsize*xsize*ysize])
	{}

	// Implementation of the destructor
	template< typename Type, uint Cellsize > 
		Grid<Type,Cellsize>::~Grid(){
		delete [] data_;
		}

	// Implementation of the function call operator
	template< typename Type, uint Cellsize >
		inline Type& 
		Grid<Type,Cellsize>::operator()( uint x, uint y, uint f )
		{
			//assert( x < xsize_ && y < ysize_ && f < Cellsize );
			return data_[y*xsize_*Cellsize+x*Cellsize+f];
		}

	// Implementation of the const function call operator// 
	template< typename Type, uint Cellsize >
		inline Type 
		Grid<Type,Cellsize>::operator()( uint x, uint y, uint f ) const
		{
			//assert( x < xsize_ && y < ysize_ && f < Cellsize );
			return data_[y*xsize_*Cellsize+x*Cellsize+f];
		}


	template< typename Type, size_t N >
	inline void swap(Grid<Type,N>& a,Grid<Type,N>& b ) /* throw() */
	{
		a.swap(b);
	}


	// Partial template specialization for Cellsize = 1
	template< typename Type >
		class Grid<Type,1>
		{
			public:
				inline Grid( uint xsize, uint ysize );
				inline ~Grid();

				inline Type& operator()( uint x, uint y);
				inline Type  operator()( uint x, uint y) const;


			private:
				uint xsize_;  // Number of nodes in x-dimension
				uint ysize_;  // Number of nodes in y-dimension
				Type* data_;  // Linearized, 1-dimensional representation of the 2D data grid
		};

	// Implementation of the initialization constructor
	template< typename Type>
		Grid<Type,1>::Grid( uint xsize, uint ysize ):
			xsize_(xsize),
			ysize_(ysize),
			data_(new Type[xsize*ysize])
	{}

	// Implementation of the destructor
	template< typename Type > 
		Grid<Type,1>::~Grid(){
		delete [] data_;
		}

	// Implementation of the function call operator
	template< typename Type>
		inline Type& 
		Grid<Type,1>::operator()( uint x, uint y)
		{
			//assert( x < xsize_ && y < ysize_ );
			return data_[y*xsize_+x];
		}

	// Implementation of the const function call operator//
	template< typename Type>
		inline Type 
		Grid<Type,1>::operator()( uint x, uint y ) const
		{
			//assert( x < xsize_ && y < ysize_ );
			return data_[y*xsize_+x];
		}


	// Partial template specialization for Cellsize = 0
	// // No class definition => compile time error
	template< typename Type >
		class Grid<Type,0>;

	typedef Grid<double,9>  PDF_Field;
	typedef Grid<double,2>  V_Field;
	typedef Grid<double,1>  D_Field;
	typedef Grid<uint,1>    Flags;
}
