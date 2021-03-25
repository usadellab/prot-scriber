/// This module contains the unit tests for struct Query. Moved the tests to a different file,
/// because they grew too large.

#[cfg(test)]
mod tests {
    use crate::default::SPLIT_DESCRIPTION_REGEX;
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
    fn test_to_similarity_matrix() {
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.qlen = 100;
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.find_bitscore_scaling_factor();
        let sim_mtrx = query.to_similarity_matrix();
        // println!("Query.to_similarity_matrix() yields:\n{:?}\n", sim_mtrx);
        let expected = arr2(&[
            [0.0, h1.similarity(&h2, &query.qlen)],
            [h1.similarity(&h2, &query.qlen), 0.0],
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
        query_frivolous.qlen = 100;
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

    #[test]
    fn test_cluster_hits() {
        // Test query that should result in a single cluster:
        let h1 = Hit::new(
            "Hit_One", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h2 = Hit::new(
            "Hit_Two", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let mut query = Query::from_qacc("Query_Curious".to_string());
        query.qlen = 100;
        query.add_hit(&h1);
        query.add_hit(&h2);
        query.cluster_hits();
        assert_eq!(query.clusters.len(), 1);
        // Assure cluster is as expected and also its sort order:
        let expected_query_cluster = vec!["Hit_One", "Hit_Two"];
        // println!(
        //     "score(h1): {}, score(h2): {}",
        //     query.hits.get("Hit_One").unwrap().query_similarity_score,
        //     query.hits.get("Hit_Two").unwrap().query_similarity_score
        // );
        assert_eq!(query.clusters.get(0).unwrap().hits, expected_query_cluster);
        assert!(
            query.hits.get("Hit_One").unwrap().query_similarity_score
                > query.hits.get("Hit_Two").unwrap().query_similarity_score
        );
        // Test query that has NO hits:
        let mut query_no_hits = Query::from_qacc("Query_So_Lonely".to_string());
        query_no_hits.qlen = 123;
        query_no_hits.cluster_hits();
        assert_eq!(query_no_hits.clusters.len(), 0);
        // Test query that should produce TWO clusters:
        let h4 = Hit::new(
            "Hit_Four", "100", "1", "50", "200", "51", "110", "500.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h5 = Hit::new(
            "Hit_Five", "100", "51", "100", "300", "201", "300", "10.0",
            "sp|P15538|C11B1_HUMAN Cytochrome P450 11B1, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11B1 PE=1 SV=5"
            );
        let mut query_frivolous = Query::from_qacc("Query_Frivolous".to_string());
        query_frivolous.qlen = 100;
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.cluster_hits();
        // println!(
        //     "Query with three Hits yields these clusters:\n{:?}",
        //     query_frivolous.clusters
        // );
        // Test includes assertion of correct sort order of clusters:
        assert_eq!(query_frivolous.clusters.len(), 2);
        let expected_query_frivolous_first_cluster = vec!["Hit_Four", "Hit_Three"];
        assert_eq!(
            query_frivolous.clusters.get(0).unwrap().hits,
            expected_query_frivolous_first_cluster
        );
        assert!(
            query_frivolous
                .hits
                .get("Hit_Four")
                .unwrap()
                .query_similarity_score
                > query_frivolous
                    .hits
                    .get("Hit_Three")
                    .unwrap()
                    .query_similarity_score
        );
        assert!(query_frivolous
            .clusters
            .get(1)
            .unwrap()
            .contains(&"Hit_Five".to_string()));
    }

    #[test]
    fn test_calculate_cluster_score() {
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "100", "51", "100", "500.0",
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
        query_frivolous.qlen = 100;
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.cluster_hits();
        // println!(
        //     "Query with three Hits yields these clusters:\n{:?}",
        //     query_frivolous.clusters
        // );
        let h3_h4_indx = query_frivolous
            .clusters
            .iter()
            .position(|clstr| clstr.contains(&"Hit_Three".to_string()))
            .unwrap();
        let cluster_h3_h4_score = query_frivolous.calculate_cluster_score(
            &query_frivolous
                .clusters
                .get(h3_h4_indx)
                .unwrap()
                .hits
                .clone(),
        );
        let expected_h3_h4_clstr_score = F64((query_frivolous.similarity(&h3) + 2.0 / 3.0) / 2.0);
        // println!(
        //     "Cluster of Hits 'Hit_Three' and 'Hit_Four' of Query '{}' has index {} and cluster-score {}",
        //     query_frivolous.id, h3_h4_indx, cluster_h3_h4_score
        // );
        assert_eq!(cluster_h3_h4_score, expected_h3_h4_clstr_score);
        let h5_indx = query_frivolous
            .clusters
            .iter()
            .position(|clstr| clstr.contains(&"Hit_Five".to_string()))
            .unwrap();
        let cluster_h5_score = query_frivolous
            .calculate_cluster_score(&query_frivolous.clusters.get(h5_indx).unwrap().hits.clone());
        let expected_h5_clstr_score = F64((query_frivolous.similarity(&h5) + 1.0 / 3.0) / 2.0);
        // println!(
        //     "Cluster of Hit 'Hit_Five' of Query '{}' has index {} and cluster-score {}",
        //     query_frivolous.id, h5_indx, cluster_h5_score
        // );
        assert_eq!(cluster_h5_score, expected_h5_clstr_score);

        // Assure that the hits have their fields `query_similarity_score` set after execution of
        // calculate_cluster_score:
        assert!(query_frivolous
            .hits
            .iter()
            // Note, that the `default` query_similarity_score is 0.0, which is why we test > 0.0:
            .all(|(_, h)| h.query_similarity_score.0 > 0.0));
    }

    #[test]
    fn test_cluster_consensus_description() {
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "50", "100", "51", "100", "500.0",
            "sp|C0LGP4|Y3475_ARATH LRR receptor-like serine/threonine-protein kinase OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
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
        query_frivolous.qlen = 100;
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        query_frivolous.cluster_hits();
        let best_scoring_cluster_consensus_description =
            query_frivolous.cluster_consensus_description(0, &(*SPLIT_DESCRIPTION_REGEX));
        // println!(
        //     "query_frivolous'\n{:?}\nbest scoring cluster consensus description is:\n{}",
        //     query_frivolous, best_scoring_cluster_consensus_description
        // );
        assert_eq!(
            best_scoring_cluster_consensus_description,
            "LRR receptor-like serine/threonine-protein kinase"
        );
    }

    #[test]
    fn test_cluster_aligned_query_region() {
        let h3 = Hit::new(
            "Hit_Three", "100", "1", "45", "100", "51", "100", "500.0",
            "sp|C0LGP4|Y3475_ARATH LRR receptor-like serine/threonine-protein kinase OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h4 = Hit::new(
            "Hit_Four", "100", "10", "50", "200", "101", "200", "100.0",
            "sp|C0LGP4|Y3475_ARATH Probable LRR receptor-like serine/threonine-protein kinase At3g47570 OS=Arabidopsis thaliana OX=3702 GN=At3g47570 PE=2 SV=1"
            );
        let h5 = Hit::new(
            "Hit_Five", "100", "51", "100", "300", "201", "300", "10.0",
            "sp|P15538|C11B1_HUMAN Cytochrome P450 11B1, mitochondrial OS=Homo sapiens OX=9606 GN=CYP11B1 PE=1 SV=5"
            );
        let mut query_frivolous = Query::from_qacc("Query_Frivolous".to_string());
        query_frivolous.qlen = 100;
        query_frivolous.add_hit(&h3);
        query_frivolous.add_hit(&h4);
        query_frivolous.add_hit(&h5);
        // println!("query_frivolous'\n{:?}", query_frivolous);
        let clstr_three_four_query_region = query_frivolous
            .cluster_aligned_query_region(&vec!["Hit_Three".to_string(), "Hit_Four".to_string()]);
        assert_eq!(clstr_three_four_query_region, Some((10u32, 45u32)));
        let clstr_three_five_query_region = query_frivolous
            .cluster_aligned_query_region(&vec!["Hit_Three".to_string(), "Hit_Five".to_string()]);
        assert_eq!(clstr_three_five_query_region, None);
    }
}
