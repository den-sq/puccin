#![allow(dead_code)]
extern crate regex;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;
extern crate indexmap;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate clap;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::fs::read_dir;
use std::path::PathBuf;
use std::path;
use std::io;

use phylo::seq::{SeqData};
use phylo::paml::PAML;
use phylo::prot::{ExtendedProteinStructure, ProteinData};
use phylo::analysis::{ProteinAnalytics, StatGroupExport, ResFeature, StatGroup};
use phylo::utils::{arithmatic_median, parse_name_list, write_serializable_data};
use phylo::clade::Clade;

use clap::{App, ArgMatches, Shell};

use serde::{Serialize};
use indexmap::IndexMap;

#[derive(Serialize, Deserialize, Debug)]
struct Conf {
	data: HashMap<String, String>
}

type COMPND = HashMap<String, Vec<String>>;
type NetCounts = HashMap<String, HashMap<String, u64>>;
type ExportStats = HashMap<String, HashMap<String, HashMap<String, StatGroupExport>>>;


fn get_data_from_filename(path: &PathBuf) -> Result<(String, String, String), io::Error> {
	let error = io::Error::new(io::ErrorKind::InvalidInput, "Invalid JSON Filename");
	match path.file_name() {
		Some(fname) => {
			match fname.to_str() {
				Some(name_str) => { 
					match name_str.split('.').next() {
						Some(file) => {
							let mut details = file.split('_');
							let pdb_number = match details.next() { Some(x) => x, None => return Err(error) };
							let taedfilenumber = match details.next() { Some(x) => x, None => return Err(error) };
							let directory = match details.next() { Some(x) => x, None => return Err(error) };
							Ok((pdb_number.to_string(), taedfilenumber.to_string(), directory.to_string()))
						},
						None => Err(error)
					}
				}, None => Err(error)
			}
		},
		None => Err(error)
	}
}

#[allow(non_snake_case)]
fn diff_paml_dNdS_data(paml_duplicates: Vec<Vec<PAML>>) -> (Vec<f64>, Vec<f64>, usize) {
	let mut gaps: Vec<f64> = Vec::new();
	let mut changes: Vec<f64> = Vec::new();
	let mut flipped: usize = 0;
	for paml_set in paml_duplicates {
		for i in 0..paml_set[0].records.len() {
			// Minimum dS Value for accurate dNdS
			//if paml_set[0].records[i].dS >= 0.0001 {
				let temp = (paml_set[0].records[i].dNdS() - paml_set[1].records[i].dNdS()).abs();
				gaps.push(temp);
				changes.push(temp / paml_set[0].records[i].dNdS());

				if paml_set[0].records[i].dNdS() < 1.0 {
					if paml_set[1].records[i].dNdS() >= 1.0 {
						flipped += 1;
						println!("!{}!{}!{}", paml_set[0].records[i].dNdS(), paml_set[1].records[i].dNdS(), paml_set[0].records[i].dS);
					}
				} else if paml_set[1].records[i].dNdS() < 1.0 {
					flipped += 1;
					println!("!{}!{}!{}", paml_set[0].records[i].dNdS(), paml_set[1].records[i].dNdS(), paml_set[0].records[i].dS);
				}

				println!("{}:{}", paml_set[0].records[i].dNdS(), temp);
			//}
		}
	}
	(gaps, changes, flipped)
}

#[allow(non_snake_case)]
fn dump_paml_dNdS_data(paml_duplicates: &[Vec<PAML>]) {
	let mut run: HashMap<&str, Vec<f64>> = 
		[	("A Clipped", Vec::new()), 
			("B Clipped", Vec::new()),
			("A", Vec::new()),
			("B", Vec::new())].iter().cloned().collect();

	for paml_set in paml_duplicates {
		for i in 0..paml_set[0].records.len() {
			if (paml_set[0].records[i].dNdS() >= 0.0000) && (paml_set[1].records[i].dNdS() >= 0.0000)  {
				if (paml_set[0].records[i].dS >= 0.0000) && (paml_set[1].records[i].dS >= 0.0000) {			
					run.entry("A Clipped").and_modify(|x| x.push(paml_set[0].records[i].dNdS()));
					run.entry("B Clipped").and_modify(|x| x.push(paml_set[1].records[i].dNdS()));
				}
				run.entry("A").and_modify(|x| x.push(paml_set[0].records[i].dNdS()));
				run.entry("B").and_modify(|x| x.push(paml_set[1].records[i].dNdS()));
			}
		}
	}

	write_serializable_data(&run, "data/raw/paml_duplicates.json");
}

fn compare_paml_runs(root_dir: String, paml_runs: Vec<(&str, &str)>) {
	let mut paml_records: Vec<Vec<PAML>> = Vec::new();

	for run in paml_runs {
		let mut paml_one = match PAML::read(&format!("{}/{}/{}.subTree_FREE.1.paml", root_dir, run.0, run.1)) {
			Ok(res) => (res),
			Err(e) => {
				eprintln!("Branch Data 1 Not Added For {}/{}:\n{}", run.0, run.1, e);
				continue
			}
		};

		paml_one.regenerate_records();

		let mut paml_two = match PAML::read(&format!("{}/{}/{}.subTree_FREE.2.paml", root_dir, run.0, run.1)) {
			Ok(res) => (res),
			Err(e) => {
				eprintln!("Branch Data 2 Not Added For {}/{}:\n{}", run.0, run.1, e);
				continue
			}
		};

		paml_two.regenerate_records();

		paml_records.push(vec![paml_one, paml_two]);
	}

	dump_paml_dNdS_data(&paml_records);
	let (mut gaps, mut changes, flipped) = diff_paml_dNdS_data(paml_records);

	println!("Median dNdS Change between PAML runs:\n {} ({}%) \n {} ({}%) Flipped", // Mean:\n {} ({}%).\n
				arithmatic_median(&mut gaps).unwrap(), arithmatic_median(&mut changes).unwrap() * 100.0, 
				flipped, (flipped as f64 / gaps.len() as f64) * 100.0); // 				arithmatic_mean::<I, f64>(gaps.into_iter()), arithmatic_mean(changes.into_iter()) * 100.0, 
	println!("Changes at 75, 90, 95 percentile:\n {}, {}, {}\n {}%, {}%, {}%",
				gaps[(gaps.len() as f64 * 0.75) as usize], gaps[(gaps.len() as f64 * 0.90) as usize], gaps[(gaps.len() as f64 * 0.95) as usize],
				changes[(changes.len() as f64 * 0.75) as usize], changes[(changes.len() as f64 * 0.90) as usize], changes[(changes.len() as f64 * 0.95) as usize]
			);

	for i in 0..gaps.len() {
		if gaps[i] > 0.0000 { println!("{}% Have No Change", (((i - 1) * 100) as f64) / gaps.len() as f64); break; }
	}
	// 4700416|56871
	println!("{:?}\n{:?}\n", gaps, changes);

}	// 4700416|56871

/// Manages the tree manipulation steps
/// 
/// Parameters:
/// 
/// # matches:  Clap argument matches.
fn manipulate_tree(matches: &ArgMatches) {
	let mut tree: Clade = Clade::read(
		matches.value_of("TREE").unwrap(), 
		matches.value_of("format").unwrap_or("extension")).unwrap();
	
	if matches.is_present("annotate") {
		tree.annotate_from_csv(matches.value_of("annotate").unwrap(), ",").unwrap(); 
	}

	if matches.is_present("trim") {
		let node_map = parse_name_list(matches.value_of("trim").unwrap()).unwrap();
		tree.trim(&node_map);
	}

	if matches.is_present("common_ancestor") {
		let mut name_list = parse_name_list(matches.value_of("common_ancestor").unwrap())
								.unwrap().keys().cloned().collect();

		tree = match tree.find_common_ancestor(&mut name_list) {
			Some(new_tree) => new_tree.clone(), // This is very expensive, need to rework
												// (Borrow checker is tricky here)
			None => {
				println!("No common ancestor found.  Species Not Identified:{:?}", name_list);
				return;
			}
		};
	}

	if matches.is_present("leaf_export") {
		if let Err(e) = tree.write_leaves_csv(matches.value_of("leaf_export").unwrap()) {
			println!("Error Writing Leaf Nodes to File: {}", e);
		}
	}

	if let Err(e) = tree.write(matches.value_of("output").unwrap(), "extension") {
		println!("Error writing tree data to file: {}", e);	
		// ""
	}
}

fn manage_bootstrap(matches: &ArgMatches, verbosity: u64) {
	let rev_format_map: HashMap<String, String> = [
		("Helix".to_string(), "A".to_string()), ("All Sites".to_string(), "*".to_string()),
		("3₁₀ Helix".to_string(), "G".to_string()), ("α-Helix".to_string(), "H".to_string()), ("π-Helix".to_string(), "I".to_string()), 
		("Turn".to_string(), "T".to_string()), ("Bend".to_string(), "S".to_string()), ("Coil".to_string(), "C".to_string()), 
		("β-Bridge".to_string(), "B".to_string()), ("β-Sheet".to_string(), "E".to_string()), ("Unidentified".to_string(), " ".to_string()), 
		("Buried".to_string(), "P".to_string()), ("Exposed".to_string(), "F".to_string())].iter().cloned().collect();

	let letter_net_counts: NetCounts;
	let old_stats: ExportStats;
	let mut adjusted_stats: HashMap<String, HashMap<String, HashMap<String, f64>>> = HashMap::new();
	adjusted_stats.insert("lpp".to_string(), 
		[("All Sites".to_string(), HashMap::new()), 
		("Substituted Sites".to_string(), HashMap::new())].iter().cloned().collect());
	adjusted_stats.insert("lpf".to_string(), 
		[("All Sites".to_string(), HashMap::new()), 
		("Substituted Sites".to_string(), HashMap::new())].iter().cloned().collect());
	adjusted_stats.insert("lnp".to_string(), 
		[("All Sites".to_string(), HashMap::new()), 
		("Substituted Sites".to_string(), HashMap::new())].iter().cloned().collect());
	adjusted_stats.insert("lnf".to_string(), 
		[("All Sites".to_string(), HashMap::new()), 
		("Substituted Sites".to_string(), HashMap::new())].iter().cloned().collect());

	// Load Stats
	match File::open(matches.value_of("stats").unwrap()) { // "stats.json"
		Ok(file) => {
			let buf_reader = BufReader::new(file);
			old_stats = match serde_json::from_reader(buf_reader) {
				Ok(analytics) => analytics,
				Err(e) => panic!("Could Not Parse Stored Pre-Calculated Data for Re-Analysis, Exiting: {}", e)
			};
		},
		Err(e) => panic!("Could Not Open Pre-Existing Pre-Calculated Data, Exiting: {}", e)
	}

	// Recaclulate Existing Stats
	for site_type in ["All Sites".to_string(), "Substituted Sites".to_string()].iter() {
		let pp_base = old_stats["Changed on a Lineage"]["Positively Selected Lineages"][site_type]["All Sites"]["P"].0;
		let pf_base = old_stats["Changed on a Lineage"]["Positively Selected Lineages"][site_type]["All Sites"]["F"].0;
		let np_base = old_stats["Changed on a Lineage"]["Negatively Selected Lineages"][site_type]["All Sites"]["P"].0;
		let nf_base = old_stats["Changed on a Lineage"]["Negatively Selected Lineages"][site_type]["All Sites"]["F"].0;

		for (key, value) in old_stats[&"Changed on a Lineage".to_string()][&"Negatively Selected Lineages".to_string()][site_type].iter() {
			adjusted_stats.get_mut("lnp").unwrap().entry(site_type.clone()).and_modify(|x| {
				x.insert(rev_format_map[key].clone(), value["P"].0 / np_base);
			});
			adjusted_stats.get_mut("lnf").unwrap().entry(site_type.clone()).and_modify(|x| {
				x.insert(rev_format_map[key].clone(), value["F"].0 / nf_base);
			});
		}
		for (key, value) in old_stats[&"Changed on a Lineage".to_string()][&"Positively Selected Lineages".to_string()][site_type].iter() {
			adjusted_stats.get_mut("lpp").unwrap().entry(site_type.clone()).and_modify(|x| {
				x.insert(rev_format_map[key].clone(), value["P"].0 / pp_base);
			});
			adjusted_stats.get_mut("lpf").unwrap().entry(site_type.clone()).and_modify(|x| {
				x.insert(rev_format_map[key].clone(), value["F"].0 / pf_base);
			});
		}
	}
	
	// load Net Counts
	match File::open(matches.value_of("counts").unwrap()) { // "dist_net_counts.json"
		Ok(file) => {
			let buf_reader = BufReader::new(file);
			letter_net_counts = match serde_json::from_reader(buf_reader) {
				Ok(analytics) => analytics,
				Err(e) => panic!("Could Not Parse Stored Pre-Calculated Data for Re-Analysis, Exiting: {}", e)
			};
		},
		Err(e) => panic!("Could Not Open Pre-Existing Pre-Calculated Data, Exiting: {}", e)
	}
	
	// Get Generated Stats
	let mut dist_net_counts: HashMap<String, IndexMap<ResFeature, u64>> = HashMap::new();

	for (base, data) in letter_net_counts.iter() {

		let mut count_vec: Vec<(&String, &u64)> = data.iter().collect();
		count_vec.sort_by(|a, b| b.1.cmp(a.1).reverse());

		dist_net_counts.insert(base.clone(), IndexMap::new());
		for (key, value) in count_vec {
			dist_net_counts.entry(base.clone())
				.and_modify(|x| { x.insert(
					ResFeature { ss: key.chars().next().unwrap(), sa: key.chars().last().unwrap() }, 
					*value as u64); }
				);
		}
	}
	let bootstrap_data: HashMap<String, Vec<StatGroup>> = 
		ProteinAnalytics::generate_bootstrap(dist_net_counts, 
			[("lp".to_string(), 56_871), ("ln".to_string(), 785_329)].iter().cloned().collect(),
			[("lp".to_string(), 4_700_416), ("ln".to_string(), 35_320_043)].iter().cloned().collect(),
			20000, vec!["lp".to_string(), "ln".to_string()]);

	// 823813|44188
	// 1145806|356568

	if verbosity > 2 { println!("{:?}", adjusted_stats); }

	// Test What Data We Have
	write_serializable_data(&adjusted_stats, "adjusted_counts.json");
	write_serializable_data(&bootstrap_data, "bootstrap_data.json");

	// Move Bootstrap Data
	let mut bootstrap_stats: HashMap<String, Vec<HashMap<String, f64>>> = 
		[("rpp".to_string(), Vec::new()),
		("rpf".to_string(), Vec::new()),
		("rnp".to_string(), Vec::new()),
		("rnf".to_string(), Vec::new())].iter().cloned().collect();

	for (base, data) in bootstrap_data {
		let bootstrap_base = format!("r{}", base.chars().last().unwrap());

		for stat_set in data {
			let mut p_data:  HashMap<String, f64> = HashMap::new();
			let mut f_data:  HashMap<String, f64> = HashMap::new();
			let mut p_total: f64 = 0.0;
			let mut f_total: f64 = 0.0;

			for (stat_name, stat) in stat_set.stats.iter() {
				let solvent_exposure: char = stat_name.chars().last().unwrap();
				if solvent_exposure == 'P' {
					p_total += stat.value;
				} else if solvent_exposure == 'F' {
					f_total += stat.value;
				}
			}
			p_total -= stat_set.stats["AP"].value;
			f_total -= stat_set.stats["AF"].value;
			if verbosity > 1 { eprintln!("{}|{}", p_total, f_total); }

			for (stat_name, stat) in stat_set.stats.iter() {
				let solvent_exposure: char = stat_name.chars().last().unwrap();

				if solvent_exposure == 'P' {
					p_data.insert(stat_name.chars().next().unwrap().to_string(),
						stat.value / p_total);
				} else if solvent_exposure == 'F' {
					f_data.insert(stat_name.chars().next().unwrap().to_string(),
						stat.value / f_total);
				}
			}
			bootstrap_stats.get_mut(&format!("{}p", bootstrap_base)).unwrap()
				.push(p_data);
			bootstrap_stats.get_mut(&format!("{}f", bootstrap_base)).unwrap()
				.push(f_data);
		}
	}
	
	let empty_results: HashMap<String, u64> =
		[("A".to_string(), 0), ("G".to_string(), 0), ("H".to_string(), 0), ("I".to_string(), 0), 
			("T".to_string(), 0), ("S".to_string(), 0), ("C".to_string(), 0), 
			("B".to_string(), 0), ("E".to_string(), 0), (" ".to_string(), 0), 
			("P".to_string(), 0), ("F".to_string(), 0)].iter().cloned().collect();

	let mut p_counts: HashMap<String, HashMap<String, u64>> = 
		[("rpp".to_string(), empty_results.clone()),
		("rpf".to_string(), empty_results.clone()),
		("rnp".to_string(), empty_results.clone()),
		("rnf".to_string(), empty_results)].iter().cloned().collect();

	for (key, results) in bootstrap_stats {
		for result in results {
			for (field, value) in result {
				let adj_str = key.replace("r", "l");
				if verbosity > 2 {
					println!("{}:{}:{}:{}:", &field, &value, adjusted_stats[&adj_str]["All Sites"][&field], adjusted_stats[&adj_str]["Substituted Sites"][&field]);
					println!("{}:{}:{}:{}:", &field,
						(value - adjusted_stats[&adj_str]["All Sites"][&field]).abs(), 
						(adjusted_stats[&adj_str]["Substituted Sites"][&field] - adjusted_stats[&adj_str]["All Sites"][&field]).abs(), 
						((value - adjusted_stats[&adj_str]["All Sites"][&field]).abs() > 
						(adjusted_stats[&adj_str]["Substituted Sites"][&field] - adjusted_stats[&adj_str]["All Sites"][&field]).abs()
						) as u64);
				}
				p_counts.get_mut(&key).unwrap().entry(field.clone()).and_modify(|x| *x +=
					((value - adjusted_stats[&adj_str]["All Sites"][&field]).abs() > 
					(adjusted_stats[&adj_str]["Substituted Sites"][&field] - adjusted_stats[&adj_str]["All Sites"][&field]).abs()
					) as u64
				);
			}
		}
	}

	if verbosity > 2 { println!("{:?}", p_counts); }

	write_serializable_data(&p_counts, "p_counts.json");

	let mut p_values: HashMap<String, HashMap<String, f64>> = 
		[("rpp".to_string(), HashMap::new()),
		("rpf".to_string(), HashMap::new()),
		("rnp".to_string(), HashMap::new()),
		("rnf".to_string(), HashMap::new())].iter().cloned().collect();

	for (key, counts) in p_counts {
		for (field, count) in counts {
			p_values.get_mut(&key).unwrap().insert(field.clone(), count as f64 / 20000_f64);
		}
	}

	if verbosity > 2 { println!("{:?}", p_values); }

	write_serializable_data(&p_values, "p_values.json");
}

fn structure_analysis(matches: &ArgMatches, verbosity: u64) {
	// Config File
	let config: Conf = match matches.value_of("config") {
		Some(config_file) => match File::open(config_file) {
			Ok(datastore) => {
				let buf_reader = BufReader::new(datastore);
				match serde_yaml::from_reader(buf_reader) {
					Ok(config_object) => config_object,
					Err(e) => panic!("Bad Config File Format: {}", e)
				}
			}, Err(e) => panic!("Config File Not Found Or Could Not Be Opened: {}", e)
		}, None => panic!("No Config File Specified.")
	};

	let mut how_far: u64 = 1;
	let mut file_list: Vec<PathBuf> = Vec::new();
	let mut analytics = ProteinAnalytics::new();
	let mut prot = ExtendedProteinStructure::new();

	let prot_lookup: COMPND = match File::open("data/interim/protein_types.json") {
		Ok(file) => {
			let buf_reader = BufReader::new(file);
			match serde_json::from_reader(buf_reader) {
				Ok(prot_lookup) => prot_lookup,
				Err(e) => panic!("Could Not Parse Stored Pre-Calculated Data for Re-Analysis, Exiting: {}", e)
			}
		},
		Err(e) => {
			panic!("Could not find protein types: {}", e)
		}
	};

	if !path::Path::new("data/interim/pre_calc_data.json").exists() {
		let dssp_file_dir = read_dir(config.data["prot"].clone()).unwrap();
		for path in dssp_file_dir {
			file_list.push(path.unwrap().path());
		}
		file_list.sort();
		let mut last_pdb = String::from("Fresh");
		for path in file_list {
			let path_string = path.to_str().unwrap();
			match File::open(path_string) {
				Ok(file) => {
					//eprintln!("At {} out of 74468", how_far);
					how_far += 1;
					let pdb_number: String; let taedfilenumber: String; let directory: String;

					match get_data_from_filename(&path) {
						Ok((a, b, c)) => {
							pdb_number = a;
							taedfilenumber = b;
							directory = c;
						},
						Err(e) => {
							eprintln!("Bad JSON Filename: {}: {}", path_string, e);
							continue;
						}
					};

					if !prot_lookup[&pdb_number].is_empty() {
						// It's a membrane protein and/or a multimeric protein
						if verbosity > 1 { println!("{} is a multimeric protein", pdb_number); }
						continue;
					}

					if directory == "2" { continue }

					if last_pdb != pdb_number {
						if last_pdb != "Fresh" {
							analytics.proteins.push(prot.clone());
						}	
						
						let buf_reader = BufReader::new(file);
						prot = match serde_json::from_reader(buf_reader) {
							Ok(protein) => ExtendedProteinStructure::extend(protein, pdb_number),
							Err(_) => {
								analytics.bad_datasources.entry("Bad JSON File".to_string())
									.and_modify(|x| { x.insert(format!("{}_{}_{}.json", pdb_number, taedfilenumber, directory)); })
									.or_insert_with(|| [taedfilenumber.to_string()].iter().cloned().collect());
								continue;
							}
						};

						if prot.missing_residue_check() { eprintln!("Residues are missing for PDB {}.", prot.id()) }
					}


					let mut paml_file = match PAML::read(&format!("{}/{}/{}.subTree_tree_1.paml_rooted",config.data["paml"],directory,taedfilenumber)) {
						Ok(res) => (res),
						Err(e) => {
							eprintln!("{}: Branch Data Not Added For {}:\n{}", how_far, taedfilenumber, e);
							analytics.bad_datasources.entry("Missing Branch Data".to_string())
								.and_modify(|x| { x.insert(taedfilenumber.to_string()); })
								.or_insert_with(|| [taedfilenumber.to_string()].iter().cloned().collect());
							continue;
						}
					};
					match paml_file.append(&format!("{}/{}/{}_1.RST",config.data["paml"],directory,taedfilenumber)) {
						Some(_) => (),
						None => { 
							eprintln!("{}:{}: RST Data Not Appended", how_far, format!("{}/{}/{}_1.RST",config.data["paml"],directory,taedfilenumber));
							analytics.bad_datasources.entry("Corrupted RST Data".to_string())
								.and_modify(|x| { x.insert(taedfilenumber.to_string()); })
								.or_insert_with(|| [taedfilenumber.to_string()].iter().cloned().collect());
							continue;
						}
					};

					paml_file.regenerate_records();
					paml_file.append_prot_structure(&prot);
					
					match paml_file.get_ancestral_alignment() {
						Some(extended_alignment) => prot.add_family(paml_file.records, extended_alignment, taedfilenumber),
						None => continue
					};

					last_pdb = prot.id().to_string();
				}
				Err(e) => {
					println!("{:?}", e);
				}
			}
		}
		analytics.proteins.push(prot);

		write_serializable_data(&analytics, "data/interim/pre_calc_data.json");
	} else {
		match File::open("data/interim/pre_calc_data.json") {
			Ok(file) => {
				let buf_reader = BufReader::new(file);
				analytics = match serde_json::from_reader(buf_reader) {
					Ok(analytics) => analytics,
					Err(e) => panic!("Could Not Parse Stored Pre-Calculated Data for Re-Analysis, Exiting: {}", e)
				};
			},
			Err(e) => panic!("Could Not Open Pre-Existing Pre-Calculated Data, Exiting: {}", e)
		}
		
	}

	analytics.dNdS_record();

	write_serializable_data(&analytics.dNdS_distribution, "data/results/dNdSvsdS.json");

	analytics.total();
	analytics.build_stats();
	analytics.generate_sample();

	write_serializable_data(&analytics.export_stats(), "data/results/stats.json");
	write_serializable_data(&analytics.export_counts(), "data/results/counts.json");
	write_serializable_data(&analytics.export_raw_counts(), "data/raw/compilation.json");
	write_serializable_data(&analytics.bad_datasources, "data/results/error_records.json");

	if verbosity > 2 { println!("{}", analytics); }
	if verbosity > 0 { 
		eprintln!("Metadata: {:?}\n{:?}", analytics.export_counts(), analytics.bad_datasources);
	}
}
	// 'B' => "ASX", 'X' => "XAA", 'Z' => "GLX", 'J' => "XLE", 'U' => "SEC", 'O' => "PYL"

fn main() {
	let yaml = load_yaml!("cli.yaml");
	let matches = App::from_yaml(yaml).get_matches();

	let verbosity: u64 = matches.occurrences_of("verbose");

	match matches.subcommand() {
        ("completions", Some(matches)) => {
            let shell = matches.value_of("SHELL").unwrap();
            App::from_yaml(yaml).gen_completions_to(
                "puccin", 
                shell.parse::<Shell>().unwrap(), 
                &mut io::stdout()
            );
		},
		("tree", Some(matches)) => { manipulate_tree(matches); },
		("bootstrap", Some(matches)) => { manage_bootstrap(matches, verbosity); },
		("paml_analysis", Some(matches)) => {
			compare_paml_runs(matches.value_of("compare_dir").unwrap().to_string(), vec![
				("157", "1245"), ("161", "4861"), ("247", "18791"), ("278", "6124"), ("289", "3757"), ("294", "14775"), ("37", "21452"), ("379", "11460"), ("394", "75300"), ("398", "12732"), ("44", "13326"), ("48", "45391"), ("491", "5029"), ("55", "20991"), ("556", "189"), ("590", "15774"), ("602", "11162"), ("608", "11700"), ("612", "2511"), ("648", "19486"), ("766", "5750"), ("84", "37973"), ("94", "24976")
			]);
		},
		("structure_analysis", Some(matches)) => { structure_analysis(matches, verbosity); },
		(&_,_) => {} // No other subcommands.
	};
}