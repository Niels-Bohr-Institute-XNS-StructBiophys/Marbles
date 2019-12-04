/*******************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

template <class T, size_t W, size_t H>
class Array2D {

  public:
    int width = W; /* width of the matrix */
    int height = H; /* height of the matrix */

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

    void add( unsigned int x, unsigned int y, T val ) {
      buffer[y*width + x] += val;
    }

    void resize( unsigned int x, unsigned int y ) {
      buffer.resize( x * y );
      width = y;
      height = x;
    }

    std::vector<T> get_buffer() {
      return buffer;
    }

    //I will need to check the type here!
    void copy_from( Array2D to_copy ) {
      buffer = to_copy.get_buffer();
    }

  private:
    std::vector<T> buffer;
};

template <class T, size_t W, size_t H, size_t D> //width, height and depth
class Array3D {

  public:
    int width        = W; /* number of rows */
    int height = H; /* number of columns */
    int depth  = D; /* thrid dimension */

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

    std::vector<T> get_buffer() {
      return buffer;
    }

    //I will need to check the type here!
    void copy_from( Array3D to_copy ) {
      buffer = to_copy.get_buffer();
    }

    void copy_from_vector( std::vector<T> to_copy, int wi, int he, int de ) {
      buffer = to_copy;
      width  = wi;
      height = he;
      depth  = de;
    }

  private:
    std::vector<T> buffer;
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> std::complex<double> pol(T r, T phi) {

    if( phi == 0 )
      return { r, 0. };
    else
      return { r * cos(phi), r * sin(phi) };
}
