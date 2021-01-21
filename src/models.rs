use super::default::FILTER_REGEXS;
use eq_float::F64;
use regex::Regex;
use std::collections::HashMap;

pub fn overlap(sstart: i64, send: i64, slen: i64, qstart: i64, qend: i64, qlen: i64) -> F64 {
    F64((((qend - qstart + 1) + (send - sstart + 1)) as f64) / ((qlen + slen) as f64))
}

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

#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Hit {
    pub id: String,
    pub slen: i64,
    pub sstart: i64,
    pub send: i64,
    pub bitscore: F64,
    pub overlap: F64,
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
            slen: slen_prsd,
            sstart: sstart_prsd,
            send: send_prsd,
            bitscore: F64(bitscore.parse().unwrap()),
            overlap: overlap(
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
}

#[derive(PartialEq, Eq, Debug, Clone, Default)]
pub struct Query {
    pub id: String,
    pub hits: HashMap<String, Hit>,
}

impl Query {
    pub fn from_qacc(id: String) -> Query {
        Query {
            id: id,
            hits: HashMap::<String, Hit>::new(),
        }
    }
    pub fn new() -> Query {
        Default::default()
    }
}

trait HasDistanceMeasure {
    fn distance(to: &Hit) -> F64;
}

impl HasDistanceMeasure for Query {
    fn distance(to: &Hit) -> F64 {
        F64(4.0)
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
    fn parse_hit_from_strs() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "123.4", "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        assert_eq!(h1.id, String::from("Hit_One"));
        assert_eq!(h1.bitscore, F64(123.4));
        assert_eq!(h1.sstart, 51);
        assert_eq!(h1.send, 110);
        assert_eq!(h1.slen, 200);
        assert_approx_eq!(h1.overlap.0, 0.37, 0.01);
        assert_eq!(
            h1.description,
            "Probable LRR receptor-like serine/threonine-protein kinase At3g47570"
        );
    }
}
