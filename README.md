# simple-elpmpp02

`ELP/MPP02` (Ephemeride Lunaire Parisienne / Lunar Solution) provides highly 
accurate positions for the Moon over a time span of several thousand years.
 This Rust library is a port of the original ELP/MPP02 implementation, 
 which can be found at  ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/ 

Several optimizations have been applied:
* SIMD (Single Instruction, Multiple Data)
* Ignoring small terms using the `tol` parameter

Like the original implementation, this library operates based on the dynamical frame of the mean ecliptic and equinox of J2000. The time scale of `t` is TT (Terrestrial Time).

## Usage

### Basic Usage

Here's a basic example of how to use the library:

```rust
use simple_elpmpp02;

fn main() {
    const J2000: f64 = 2451545.0; // 2000-01-01, 12:00:00 TT
    let coords = simple_elpmpp02::cartesian(J2000, 0.0);
    
    println!("Moon's position on 2000-01-01:");
    println!("X : {} km", coords.0);
    println!("Y : {} km", coords.1);
    println!("Z : {} km", coords.2);
    println!("X': {} km/d", coords.3);
    println!("Y': {} km/d", coords.4);
    println!("Z': {} km/d", coords.5);
}
```

### Scripts

* `elpmpp02.py`: Supports the following features:
  * Downloads the official ELP/MPP02 implementation and coefficients from the source.
  * Replicates the functionality of the original implementation.
  * Generates `.bin` files. You may need to modify this file if you want to reduce the bin file size.
* `generate_testcases.py`: Reads `ELPMPP02.PY.TXT` (the output of `elpmpp02.py`) and prints the test function.


## Acknowledgements

This implementation is based on the ELP/MPP02 theory developed by J. Chapront and G. Francou at the Bureau des Longitudes, Paris.