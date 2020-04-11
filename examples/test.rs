fn main() {
    let mut v: Vec<fn(usize)->usize> = Vec::new();
    v.push(|x: usize| x+1);
    println!("{}", v[0](1));
}