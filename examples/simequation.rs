use ndarray::*;
use ndarray_linalg::*;

fn main() {
    let a = arr2(&[[  2.0, -4.0,  6.0],
                   [ -1.0,  7.0, -8.0],
                   [  1.0,  1.0, -2.0]]);
    let b = arr1(&[  1.0,  0.0,  3.0]);
    println!("*** simultaneous equations ***");
    println!(" 2x - 4y + 6z =  1");
    println!("- x + 7y - 8z =  0");
    println!("  x +  y - 2z =  3");
    println!("** solution **");
    let x = a.solve(&b).unwrap();
    println!("x = {:.6}", x[0]);
    println!("y = {:.6}", x[1]);
    println!("z = {:.6}", x[2]);
    println!("");

    let a = arr2(&[[  1.0,  1.0, -1.0],
                   [  5.0,  2.0, -6.0],
                   [  4.0,  1.0, -5.0]]);
    let b = arr1(&[  2.0,  5.0,  3.0]);
    println!("*** simultaneous equations ***");
    println!("  x +  y -  z =  2");
    println!(" 5x + 2y - 6z =  5");
    println!(" 4x +  y - 5z =  3");
    println!("** solution **");
    let x = a.solve(&b).unwrap();
    // [x, y, z] = 1/3*[1, 5, 0] + s*[4, -1, 3]
    println!("x = {:.6}", x[0]);
    println!("y = {:.6}", x[1]);
    println!("z = {:.6}", x[2]);
    println!("");
    
    let a = arr2(&[[  1.0,  2.0, -2.0],
                   [  2.0, -2.0,  2.0],
                   [  1.0,  6.0, -6.0]]);
    let b = arr1(&[  2.0,  3.0, -1.0]);
    println!("*** simultaneous equations ***");
    println!("  x + 2y - 2z =  2");
    println!(" 2x - 2y + 2z =  3");
    println!("  x + 6y - 6z = -1");
    println!("** solution **");
    let x = a.solve(&b);
    // no solution
    println!("{:?}", x);
    println!("");

    let a = arr2(&[[  3.0, -6.0,  3.0,  9.0],
                   [  1.0,  0.0,  7.0, -5.0],
                   [ -1.0,  6.0,  7.0, -9.0],
                   [  0.0,  2.0,  2.0,  6.0]]);
    let b = arr1(&[6.0, 6.0, -4.0, -10.0]);
    println!("*** simultaneous equations ***");
    println!(" 3x - 6y + 3z + 9u =   6");
    println!("  x      + 7z - 5u =   6");
    println!("- x + 6y + 7z - 9u = - 4");
    println!("      2y + 2z + 6u = -10");
    let f = a.factorize_into().unwrap();
    let x = f.solve_into(b).unwrap();
    println!("*** solve by LU decomposition ***");
    println!("x = {:.6}", x[0]);
    println!("y = {:.6}", x[1]);
    println!("z = {:.6}", x[2]);
    println!("u = {:.6}", x[3]);
}