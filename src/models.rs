use super::default::FILTER_REGEXS;
use super::default::SPLIT_DESCRIPTION_REGEX;
use eq_float::F64;
use regex::Regex;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::collections::HashSet;

pub fn filter_stitle(stitle: &str, regexs: &Vec<Regex>) -> String {
    String::from(
        regexs
            .iter()
            .fold(stitle.to_string(), |accumulated, current| {
                current.replace_all(&accumulated, "").to_string()
            })
            .trim(),
    )
}

pub fn overlap_with_query(
    sstart: i64,
    send: i64,
    slen: i64,
    qstart: i64,
    qend: i64,
    qlen: i64,
) -> F64 {
    F64((((qend - qstart + 1) + (send - sstart + 1)) as f64) / ((qlen + slen) as f64))
}

#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Hit {
    pub id: String,
    pub qstart: i64,
    pub qend: i64,
    pub bitscore: F64,
    pub overlap_with_query: F64,
    pub description: String,
}

impl Hit {
    pub fn new(
        id: &str,
        qlen: &str,
        qstart: &str,
        qend: &str,
        slen: &str,
        sstart: &str,
        send: &str,
        bitscore: &str,
        stitle: &str,
    ) -> Hit {
        let slen_prsd = slen.parse().unwrap();
        let sstart_prsd = sstart.parse().unwrap();
        let send_prsd = send.parse().unwrap();
        let qlen_prsd = qlen.parse().unwrap();
        let qstart_prsd = qstart.parse().unwrap();
        let qend_prsd = qend.parse().unwrap();
        Hit {
            id: String::from(id),
            qstart: qstart_prsd,
            qend: qend_prsd,
            bitscore: F64(bitscore.parse().unwrap()),
            overlap_with_query: overlap_with_query(
                sstart_prsd,
                send_prsd,
                slen_prsd,
                qstart_prsd,
                qend_prsd,
                qlen_prsd,
            ),
            description: filter_stitle(&stitle, &(*FILTER_REGEXS)),
        }
    }

    /// Overlap between hits `o.hits( h1, h2 )` is defined as the overlap between the query's
    /// region the two respective hits aligned to. _Note_: Make sure both hits belong to the same
    /// query! Returns the overlap (`f64`) calculated as follows:
    ///
    /// o.hits( h1, h2 ) =
    ///   ( min(qend.h1, qend.h2) - max(qstart.h1, qstart.h2) ) / qlen
    ///
    /// _Note_ that the above `qlen` is of course identical between the two hits (`self` and
    /// `with`).
    ///
    /// # Arguments
    ///
    /// * `&self` - reference to an instance of Hit
    /// * `with: &Hit` - the hit which which to calculate the overlap
    /// * `qlenL &f64` - the query sequence's length
    fn overlap_with_hit(&self, with: &Hit, qlen: &f64) -> f64 {
        (min(self.qend, with.qend) as f64 - max(self.qstart, with.qstart) as f64) / qlen
    }

    fn description_words(&self) -> HashSet<&str> {
        (*SPLIT_DESCRIPTION_REGEX)
            .split(&self.description)
            .into_iter()
            .collect::<HashSet<&str>>()
    }

    fn description_distance(&self, to: &Hit) -> f64 {
        let self_words = self.description_words();
        let to_words = to.description_words();
        let intersection: HashSet<_> = self_words.intersection(&to_words).collect();
        intersection.len() as f64 / min(self_words.len(), to_words.len()) as f64
    }

    fn distance(&self, to: &Hit, qlen: &f64) -> f64 {
        let o = self.overlap_with_hit(to, qlen);
        let dd = self.description_distance(to);
        1.0 - (o + dd) / 2.0
    }
}

#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Query {
    pub id: String,
    pub qlen: F64,
    pub bitscore_scaling_factor: F64,
    pub hits: HashMap<String, Hit>,
}

impl Query {
    /// Creates an empty Query with an initialized empty HashMap (`hits`) and initialized `Default`
    /// `bitscore_scaling_factor`. Sets the `id` field to the provided argument.
    ///
    /// # Arguments
    ///
    /// * id - The `qacc` value, i.e. the sequence accession (identifier) of the query sequence
    /// as provided to the used sequence similarity search tool (e.g. Blast or Diamond).
    pub fn from_qacc(id: String) -> Query {
        Query {
            id: id,
            qlen: Default::default(),
            bitscore_scaling_factor: Default::default(),
            hits: HashMap::<String, Hit>::new(),
        }
    }

    /// Creates and empty (`Default`) instance of struct Query.
    pub fn new() -> Query {
        Default::default()
    }

    /// Adds a hit instance to the query's (`&self`) hits HashMap field. The hit is added only, if
    /// it either does not yet exist among `hits` or if its biscore is higher than the recorded hit
    /// of identical `id`. In the latter case the argument hit replaces the previously stored one.
    /// Function returns the number of hits in the query's `hits` HashMap.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - mutable reference to self (the query)
    /// * `hit` - a reference to the hit instance to be added
    pub fn add_hit(&mut self, hit: &Hit) -> usize {
        if !self.hits.contains_key(&hit.id)
            || self.hits.get(&hit.id).unwrap().bitscore < hit.bitscore
        {
            self.hits.insert(hit.id.clone(), hit.clone());
        }
        self.hits.len()
    }

    /// Identifies the bitscore-scaling-factor as the maximum of all bitscores of the query
    /// (`&self`) hits.
    ///
    /// # Arguments
    ///
    /// * `&mut self` - mutable reference to self (the query)
    pub fn find_bitscore_scaling_factor(&mut self) {
        // find maximum bitscore
        self.bitscore_scaling_factor = F64(self.hits.iter().fold(0.0, |accumulated, (_, hit)| {
            if hit.bitscore.0 > accumulated {
                hit.bitscore.0
            } else {
                accumulated
            }
        }));
    }

    /// Calculates the distance between the query (`self`) and the argument hit (`to`) as the mean
    /// of overlap_with_query and bitscore, where the bitscore is scaled by dividing it with the max. bitscore
    /// found for all the query's hits. Returns the distance ranging from zero to one, inclusive.
    ///
    /// # Arguments
    ///
    /// * `&self` - the query
    /// * `to: &Hit` - the hit (one of the query's) to which the calculate the distance to
    fn distance(&self, to: &Hit) -> f64 {
        if !self.hits.contains_key(&to.id) {
            panic!(
                "The query '{}' does not contain the hit '{}' to which calculate the distance to.",
                self.id, to.id
            );
        }
        1.0 - (to.overlap_with_query.0 + to.bitscore.0 / self.bitscore_scaling_factor.0) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn default_filter_regexs_extract_uni_prot_descriptions() {
        let t1 = "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1";
        assert_eq!(
            filter_stitle(t1, &(*FILTER_REGEXS)),
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }

    #[test]
    fn query_add_hit_only_uses_higest_bitscore() {
        let high = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let low = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "1.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let highest = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "666.6",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let other = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&high);
        assert_eq!(query.hits.len(), 1);
        query.add_hit(&low);
        assert_eq!(query.hits.len(), 1);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), high);
        query.add_hit(&other);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), high);
        assert_eq!(*query.hits.get("Hit_Two").unwrap(), other);
        query.add_hit(&highest);
        assert_eq!(query.hits.len(), 2);
        assert_eq!(*query.hits.get("Hit_One").unwrap(), highest);
        assert_eq!(*query.hits.get("Hit_Two").unwrap(), other);
    }

    #[test]
    fn test_calculate_bitscore() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "666.6",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "51", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        assert_eq!(query.bitscore_scaling_factor.0, 666.6);
    }

    #[test]
    fn parse_hit_from_strs() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(h1.id, String::from("Hit_One"));
        assert_eq!(h1.bitscore, F64(123.4));
        assert_eq!(h1.qstart, 1);
        assert_eq!(h1.qend, 50);
        assert_approx_eq!(h1.overlap_with_query.0, 0.37, 0.01);
        assert_eq!(
            h1.description,
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }

    #[test]
    fn calculate_overlap_between_hits() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "25", "50", "200", "76", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let o = h1.overlap_with_hit(&h2, &100.0);
        assert_eq!(o, 0.25);
    }

    #[test]
    fn split_description_into_words() {
        let h = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let expected_words: HashSet<&str> = vec![
            "Probable",
            "LRR",
            "receptor-like",
            "serine/threonine-protein",
            "kinase",
            "At3g47570",
        ]
        .into_iter()
        .collect();
        assert_eq!(h.description_words(), expected_words);
    }

    #[test]
    fn distance_between_hits() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "25", "50", "200", "76", "110", "123.4",
            "sp|C0LGP4|Y3475_ARATH LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let d = h1.distance(&h2, &100.0);
        assert_eq!(d, 0.375);
    }

    #[test]
    fn distance_between_query_and_hit() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        assert_eq!(h2.overlap_with_query.0, 0.5);
        assert_eq!(query.bitscore_scaling_factor.0, 500.0);
        let d = query.distance(&h2);
        // distance( query, h2 ) =
        // 1.0 - likeliness-score =
        // 1.0 - mean( overlap[ .5 ] + scaled_bitscore[ 1/5 ] )
        assert_eq!(d, (1.0 - (0.5 + 1.0 / 5.0) / 2.0));
    }
}
