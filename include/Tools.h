template <class T, size_t W, size_t H>
class Array2D {

  public:
    const int width = W; /* width of the matrix */
    const int height = H; /* height of the matrix */

    Array2D() : buffer( width * height ) {}

    inline T at( unsigned int x, unsigned int y ) {
        return buffer[y*width + x];
    }

    void set( unsigned int x, unsigned int y, float val ) {
      buffer[y*width + x] = val;
    }

    void initialize( T value ) {
      fill( buffer.begin(), buffer.end(), value );
    }

  private:
    std::vector<T> buffer;
};

template <class T, size_t W, size_t H, size_t D> //width, height and depth
class Array3D {

  public:
    int width        = W; /* number of rows */
    const int height = H; /* number of columns */
    const int depth  = D; /* thrid dimension */

    Array3D() : buffer( width * height * depth ) {}

    inline T at( unsigned int x, unsigned int y, unsigned int z ) {
        return buffer[y*width + z*width*depth + x ];
    }

    void set( unsigned int x, unsigned int y, unsigned int z, T val ) {
      buffer[y*width + z*width*depth + x] = val;
    }

    void add( unsigned int x, unsigned int y, unsigned int z, T val ) {
      buffer[y*width + z*width*depth + x] += val;
    }

    void resize_width( int new_W ) {
      buffer.resize( depth * height * new_W );
      width = new_W;
    }

    void initialize( T value ) {
      fill( buffer.begin(), buffer.end(), value );
    }

  private:
    std::vector<T> buffer;
};
