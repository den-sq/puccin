//build.rs
#[macro_use]
extern crate clap;
 
use clap::{Shell, App};
use std::env;


fn main() {
    let yaml = load_yaml!("src/cli.yaml");
    let mut app = App::from_yaml(yaml);
    
    app.gen_completions(
        "puccin",                       // We specify the bin name manually
        Shell::Bash,                    // Which shell to build completions for
        env::var("OUT_DIR").unwrap());  // Where write the completions to
}