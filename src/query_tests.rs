/// This module contains the unit tests for struct Query. Moved the tests to a different file,
/// because they grew too large.

#[cfg(test)]
mod tests {
    use crate::hit::*;
    use crate::query::*;
    use eq_float::F64;
    use ndarray::arr2;

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
    fn similarity_between_query_and_hit() {
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
        let d = query.similarity(&h2);
        // similarity( query, h2 ) =
        // mean( overlap[ .5 ] + scaled_bitscore[ 1/5 ] )
        assert_eq!(d, (0.5 + 1.0 / 5.0) / 2.0);
    }

    #[test]
    pub fn test_to_similarity_matrix() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.qlen = F64(100.0);
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        let sim_mtrx = query.to_similarity_matrix();
        // println!("Query.to_similarity_matrix() yields:\n{:?}\n", sim_mtrx);
        let expected = arr2(&[
            [0.0, h1.similarity(&h2, &query.qlen.0)],
            [h1.similarity(&h2, &query.qlen.0), 0.0],
        ]);
        assert_eq!(sim_mtrx.1, expected);
        assert!(query.seq_sim_mtrx_node_names.len() == 2);
        assert_eq!(query.seq_sim_mtrx_node_names, vec!["Hit_One", "Hit_Two"]);
        // Test query that has hits aligning to disjoint regions of the query sequence and do not
        // share words in their respective descriptions:
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h4 = Hit::new(
            "Hit_Four", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h5 = Hit::new(
            "Hit_Five", "100", "51", "100", "300", "201", "300", "10.0",
            "sp|P15538|C11B1_HUMAN Cytochrome P450 11B1, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11B1 PE=1 SV=5"
            );
        let mut query_frivolous = Query::from_qacc("Query_Frivolous".to_string());
        query_frivolous.qlen = F64(100.0);
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        let sim_mtrx = query_frivolous.to_similarity_matrix();
        // println!(
        //     "Query with hits that align to disjoint regions and do not share words in their respective descriptions yield similarity matrix:\n{:?}",
        //     sim_mtrx
        // );
        let expected_sim_mtrx = arr2(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.75], [0.0, 0.75, 0.0]]);
        assert_eq!(sim_mtrx.1, expected_sim_mtrx);
    }
}
