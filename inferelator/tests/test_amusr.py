import unittest
import os
import pandas as pd
import numpy as np
from inferelator.regression import amusr_regression
from inferelator import amusr_workflow
from inferelator.crossvalidation_workflow import create_puppet_workflow
import numpy.testing as npt
import pandas.util.testing as pdt


class TestAMuSRWorkflow(unittest.TestCase):

    def setUp(self):
        self.expr = pd.DataFrame([[2, 28, 0, 16, 1, 3], [6, 21, 0, 3, 1, 3], [4, 39, 0, 17, 1, 3],
                                  [8, 34, 0, 7, 1, 3], [6, 26, 0, 3, 1, 3], [1, 31, 0, 1, 1, 4],
                                  [3, 27, 0, 5, 1, 4], [8, 34, 0, 9, 1, 3], [1, 22, 0, 3, 1, 4],
                                  [9, 33, 0, 17, 1, 2]],
                                 columns=["gene1", "gene2", "gene3", "gene4", "gene5", "gene6"]).transpose()
        self.meta = pd.DataFrame({"Condition": ["A", "B", "C", "C", "B", "B", "A", "C", "B", "C"],
                                  "Genotype": ['WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT']})
        self.prior = pd.DataFrame([[0, 1], [0, 1], [1, 0], [0, 0]], index=["gene1", "gene2", "gene4", "gene5"],
                                  columns=["gene3", "gene6"])
        self.gold_standard = self.prior.copy()
        self.gene_list = pd.DataFrame({"SystematicName": ["gene1", "gene2", "gene3", "gene4", "gene7", "gene6"]})
        self.tf_names = ["gene3", "gene6"]
        self.workflow = create_puppet_workflow(base_class=amusr_workflow.MultitaskLearningWorkflow,
                                               regression_class=amusr_regression.AMUSRRegressionWorkflow,
                                               result_processor_class=amusr_workflow.ResultsProcessorMultiTask)
        self.workflow = self.workflow(self.expr, self.meta, self.prior, self.gold_standard)
        self.workflow.tf_names = self.tf_names
        self.workflow.gene_metadata = self.gene_list
        self.workflow.gene_list_index = "SystematicName"
        self.workflow.create_output_dir = lambda *x: None
        self.workflow.input_dir = os.path.join(os.path.dirname(__file__), "../../data/dream4")

    def test_task_separation(self):
        self.workflow.process_priors_and_gold_standard()
        self.workflow.separate_tasks_by_metadata()
        self.assertEqual(self.workflow.n_tasks, 3)
        self.assertEqual(self.workflow.tasks_names, ["A", "B", "C"])
        self.assertEqual(list(map(lambda x: x.shape, self.workflow.expression_matrix)), [(5, 2), (5, 4), (5, 4)])
        self.assertEqual(list(map(lambda x: x.shape, self.workflow.meta_data)), [(2, 2), (4, 2), (4, 2)])

    def test_task_processing(self):
        self.workflow.startup_finish()
        self.assertEqual(self.workflow.regulators.tolist(), ["gene3", "gene6"])
        self.assertEqual(self.workflow.targets.tolist(), ["gene1", "gene2", "gene4", "gene6"])
        self.assertEqual(len(self.workflow.task_design), 3)
        self.assertEqual(len(self.workflow.task_response), 3)
        self.assertEqual(len(self.workflow.task_meta_data), 3)
        self.assertEqual(len(self.workflow.task_bootstraps), 3)
        pdt.assert_frame_equal(self.workflow.task_design[0],
                               pd.DataFrame([[16., 5.], [15., 15.]], index=["gene3", "gene6"], columns=[0, 6]),
                               check_dtype=False)
        pdt.assert_frame_equal(self.workflow.task_response[0],
                               pd.DataFrame([[2, 3], [28, 27], [16, 5], [3, 4]],
                                            index=["gene1", "gene2", "gene4", "gene6"], columns=[0, 6]),
                               check_dtype=False)

    def test_result_processor_random(self):
        beta1 = pd.DataFrame(np.array([[1, 0], [0.5, 0], [0, 1], [0.5, 0]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        beta2 = pd.DataFrame(np.array([[0, 0], [0.5, 0], [1, 0], [0.5, 0]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        beta3 = pd.DataFrame(np.array([[0.5, 0.2], [0.5, 0.1], [0.5, 0.2], [0.5, 0.2]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb1 = pd.DataFrame(np.array([[0.75, 0], [0.25, 0], [0.75, 0], [0.25, 0]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb2 = pd.DataFrame(np.array([[0, 0], [1, 0], [0, 0], [1, 0]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb3 = pd.DataFrame(np.array([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        self.workflow.startup_finish()
        self.workflow.gold_standard_filter_method = 'overlap'
        self.workflow.emit_results([[beta1, beta1], [beta2, beta2], [beta3, beta3]],
                                   [[rb1, rb1], [rb2, rb2], [rb3, rb3]],
                                   self.workflow.gold_standard,
                                   self.workflow.priors_data)
        self.assertAlmostEqual(self.workflow.aupr, 0.3416666666666667)
        self.assertEqual(self.workflow.n_interact, 0)
        self.assertEqual(self.workflow.precision_interact, 0)

    def test_result_processor_perfect(self):
        beta1 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [0.5, 0]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        beta2 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [0.5, 0]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        beta3 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [0.5, 0.2]]),
                             ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb1 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [0.25, 0]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb2 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [1, 0]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        rb3 = pd.DataFrame(np.array([[0, 1], [0, 1], [1, 0], [0.5, 0.5]]),
                           ["gene1", "gene2", "gene4", "gene6"], ["gene3", "gene6"])
        self.workflow.startup_finish()
        self.workflow.emit_results([[beta1, beta1], [beta2, beta2], [beta3, beta3]],
                                   [[rb1, rb1], [rb2, rb2], [rb3, rb3]],
                                   self.workflow.gold_standard,
                                   self.workflow.priors_data)
        self.assertAlmostEqual(self.workflow.aupr, 1)

    def test_load_data_and_split(self):
        self.workflow.read_expression(file="expression.tsv")
        self.workflow.read_metadata(file="meta_data.tsv")
        self.workflow.meta_data['Group'] = "A"
        self.workflow.meta_data.loc[[False] * 200 + [True] * 221, "Group"] = "B"
        self.workflow.meta_data_task_column = "Group"
        self.workflow.prepare_tasks()

        self.assertEqual(self.workflow.n_tasks, 2)
        self.assertListEqual(self.workflow.tasks_names, ["A", "B"])
        self.assertEqual(self.workflow.expression_matrix[0].shape, (100, 200))
        self.assertEqual(self.workflow.expression_matrix[1].shape, (100, 221))
        self.assertEqual(self.workflow.meta_data[0].shape, (200, 6))
        self.assertEqual(self.workflow.meta_data[1].shape, (221, 6))

    def test_load_taskwise_files(self):
        self.workflow.read_expression(file=["expression.tsv", "expression.tsv"])
        self.workflow.read_metadata(file=["meta_data.tsv", "meta_data.tsv"])
        self.workflow.prepare_tasks()
        self.workflow.check_tasks()

        self.assertEqual(self.workflow.n_tasks, 2)
        self.assertListEqual(self.workflow.tasks_names, ["0", "1"])
        self.assertEqual(self.workflow.expression_matrix[0].shape, (100, 421))
        self.assertEqual(self.workflow.expression_matrix[1].shape, (100, 421))
        self.assertEqual(self.workflow.meta_data[0].shape, (421, 5))
        self.assertEqual(self.workflow.meta_data[1].shape, (421, 5))

    def test_taskwise_meta_with_defaults(self):
        self.workflow.read_expression(file=["expression.tsv", "expression.tsv"])
        self.workflow.read_metadata(file=["meta_data.tsv", None])
        self.workflow.prepare_tasks()
        self.workflow.check_tasks()

        self.assertEqual(self.workflow.n_tasks, 2)
        self.assertListEqual(self.workflow.tasks_names, ["0", "1"])
        self.assertEqual(self.workflow.expression_matrix[0].shape, (100, 421))
        self.assertEqual(self.workflow.expression_matrix[1].shape, (100, 421))
        self.assertEqual(self.workflow.meta_data[0].shape, (421, 5))
        self.assertEqual(self.workflow.meta_data[1].shape, (421, 4))

    def test_extract_taskwise_metadata(self):
        mixed_data = pd.concat((self.expr.transpose(), self.meta), axis=1)
        self.workflow.expression_matrix = [mixed_data.copy(), mixed_data.copy()]
        self.workflow.n_tasks = 2
        self.workflow.expression_matrix_columns_are_genes = True
        self.workflow.extract_metadata_from_expression_matrix = True
        self.workflow.expression_matrix_metadata = ["Condition", "Genotype"]

        self.workflow.read_metadata()
        self.workflow.transpose_expression_matrix()
        self.workflow.prepare_tasks()
        self.workflow.check_tasks()

        self.assertEqual(self.workflow.n_tasks, 2)
        self.assertListEqual(self.workflow.tasks_names, ["0", "1"])
        self.assertEqual(self.workflow.expression_matrix[0].shape, (6, 10))
        self.assertEqual(self.workflow.expression_matrix[1].shape, (6, 10))
        self.assertEqual(self.workflow.meta_data[0].shape, (10, 2))
        self.assertEqual(self.workflow.meta_data[1].shape, (10, 2))

    def test_taskwise_prior_gs(self):
        self.workflow.read_priors(priors_file=["gold_standard.tsv", "gold_standard.tsv"],
                                  gold_standard_file=["gold_standard.tsv", "gold_standard.tsv"])

        self.assertTrue(isinstance(self.workflow.priors_data, list))
        self.assertTrue(isinstance(self.workflow.gold_standard, list))
        pdt.assert_frame_equal(self.workflow.priors_data[0], self.workflow.priors_data[1])
        pdt.assert_frame_equal(self.workflow.gold_standard[0], self.workflow.gold_standard[1])

    def test_taskwise_default_metadata(self):
        self.workflow.meta_data_handlers = ["nonbranching", "branching"]
        self.workflow.read_expression(file=["expression.tsv", "expression.tsv"])
        self.workflow.read_metadata(file=[None, None])
        self.assertEqual(self.workflow.meta_data[0].shape, (421, 3))
        self.assertEqual(self.workflow.meta_data[1].shape, (421, 4))


    def test_expr_meta_mismatch(self):
        self.workflow.read_expression(file="expression.tsv")

        with self.assertRaises(AssertionError):
            self.workflow.read_metadata(file=["meta_data.tsv", "meta_data.tsv"])

        self.workflow.n_tasks = 2
        self.workflow.read_metadata(file=["meta_data.tsv", "meta_data.tsv"])
        self.workflow.n_tasks = None

        with self.assertRaises(NotImplementedError):
            self.workflow.prepare_tasks()

        self.workflow.expression_matrix = [self.workflow.expression_matrix, self.workflow.expression_matrix]

        with self.assertRaises(ValueError):
            self.workflow.check_tasks()

        self.workflow.n_tasks = 2
        self.assertTrue(self.workflow.check_tasks())

        self.workflow.n_tasks = 4
        with self.assertRaises(ValueError):
            self.workflow.check_tasks()

        self.workflow.expression_matrix = self.workflow.expression_matrix + self.workflow.expression_matrix
        with self.assertRaises(ValueError):
            self.workflow.check_tasks()


class TestAMuSRrunner(unittest.TestCase):

    def test_format_priors_noweight(self):
        tfs = ['tf1', 'tf2']
        priors = [pd.DataFrame([[0, 1], [1, 0]], index=['gene1', 'gene2'], columns=tfs),
                  pd.DataFrame([[0, 0], [1, 0]], index=['gene1', 'gene2'], columns=tfs)]
        gene1_prior = amusr_regression._format_prior(priors, 'gene1', [0, 1], 1)
        gene2_prior = amusr_regression._format_prior(priors, 'gene2', [0, 1], 1)
        npt.assert_almost_equal(gene1_prior, np.array([[1., 1.], [1., 1.]]))
        npt.assert_almost_equal(gene2_prior, np.array([[1., 1.], [1., 1.]]))

    def test_format_priors_pweight(self):
        tfs = ['tf1', 'tf2']
        priors = [pd.DataFrame([[0, 1], [1, 0]], index=['gene1', 'gene2'], columns=tfs),
                  pd.DataFrame([[0, 0], [1, 0]], index=['gene1', 'gene2'], columns=tfs)]
        gene1_prior = amusr_regression._format_prior(priors, 'gene1', [0, 1], 1.2)
        gene2_prior = amusr_regression._format_prior(priors, 'gene2', [0, 1], 1.2)
        npt.assert_almost_equal(gene1_prior, np.array([[1.09090909, 1.], [0.90909091, 1.]]))
        npt.assert_almost_equal(gene2_prior, np.array([[0.90909091, 0.90909091], [1.09090909, 1.09090909]]))

    def test_sum_squared_errors(self):
        X = [np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]),
             np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])]
        Y = [np.array([3, 3, 3]),
             np.array([3, 3, 3])]
        W = np.array([[1, 0], [1, 0], [1, 0]])
        self.assertEqual(amusr_regression.sum_squared_errors(X, Y, W, 0), 0)
        self.assertEqual(amusr_regression.sum_squared_errors(X, Y, W, 1), 27)

    def test_amusr_regression(self):
        des = [np.array([[1, 1, 3], [0, 0, 2], [0, 0, 1]]).astype(float),
               np.array([[1, 1, 3], [0, 0, 2], [0, 0, 1]]).astype(float)]
        res = [np.array([1, 2, 3]).reshape(-1, 1).astype(float),
               np.array([1, 2, 3]).reshape(-1, 1).astype(float)]
        tfs = ['tf1', 'tf2', 'tf3']
        targets = ['gene1', 'gene2']
        priors = [pd.DataFrame([[0, 1, 1], [1, 0, 1]], index=targets, columns=tfs),
                  pd.DataFrame([[0, 0, 1], [1, 0, 1]], index=targets, columns=tfs)]
        runner = amusr_regression.AMuSR_regression([pd.DataFrame(des[0], columns=tfs)],
                                                   [pd.DataFrame(res[0], columns=["gene1"])],
                                                   None)
        gene1_prior = amusr_regression._format_prior(priors, 'gene1', [0, 1], 1.)
        gene2_prior = amusr_regression._format_prior(priors, 'gene2', [0, 1], 1.)
        output = []
        output.append(
            amusr_regression.run_regression_EBIC(des, res, ['tf1', 'tf2', 'tf3'], [0, 1], 'gene1', gene1_prior))
        output.append(
            amusr_regression.run_regression_EBIC(des, res, ['tf1', 'tf2', 'tf3'], [0, 1], 'gene2', gene2_prior))
        out0 = pd.DataFrame([['tf3', 'gene1', -1, 1],
                             ['tf3', 'gene1', -1, 1]],
                            index=pd.MultiIndex(levels=[[0, 1], [0]],
                                                codes=[[0, 1], [0, 0]]),
                            columns=['regulator', 'target', 'weights', 'resc_weights'])
        out1 = pd.DataFrame([['tf3', 'gene2', -1, 1],
                             ['tf3', 'gene2', -1, 1]],
                            index=pd.MultiIndex(levels=[[0, 1], [0]],
                                                codes=[[0, 1], [0, 0]]),
                            columns=['regulator', 'target', 'weights', 'resc_weights'])
        pdt.assert_frame_equal(pd.concat(output[0]), out0, check_dtype=False)
        pdt.assert_frame_equal(pd.concat(output[1]), out1, check_dtype=False)

    def test_unaligned_regression_genes(self):
        tfs = ['tf1', 'tf2', 'tf3']
        targets = ['gene1', 'gene2', 'gene3']
        targets1 = ['gene1', 'gene2']
        targets2 = ['gene1', 'gene3']
        des = [pd.DataFrame(np.array([[1, 1, 3], [0, 0, 2], [0, 0, 1]]).astype(float), columns=tfs),
               pd.DataFrame(np.array([[1, 1, 3], [0, 0, 2], [0, 0, 1]]).astype(float), columns=tfs)]
        res = [pd.DataFrame(np.array([[1, 1], [2, 2], [3, 3]]).astype(float), columns=targets1),
               pd.DataFrame(np.array([[1, 1], [2, 2], [3, 3]]).astype(float), columns=targets2)]
        priors = pd.DataFrame([[0, 1, 1], [1, 0, 1], [1, 0, 1]], index=targets, columns=tfs)

        r = amusr_regression.AMuSR_regression(des, res, tfs=tfs, genes=targets, priors=priors)

        out = [pd.DataFrame([['tf3', 'gene1', -1, 1], ['tf3', 'gene1', -1, 1]],
                            index=pd.MultiIndex(levels=[[0, 1], [0]], codes=[[0, 1], [0, 0]]),
                            columns=['regulator', 'target', 'weights', 'resc_weights']),
               pd.DataFrame([['tf3', 'gene2', -1, 1]],
                            index=pd.MultiIndex(levels=[[0, 1], [0]], codes=[[0], [0]]),
                            columns=['regulator', 'target', 'weights', 'resc_weights']),
               pd.DataFrame([['tf3', 'gene3', -1, 1]],
                            index=pd.MultiIndex(levels=[[0, 1], [0]], codes=[[1], [0]]),
                            columns=['regulator', 'target', 'weights', 'resc_weights'])]

        regress_data = r.regress()
        for i in range(len(targets)):
            pdt.assert_frame_equal(pd.concat(regress_data[i]), out[i], check_dtype=False)

    def test_process_output(self):
        cols = ["regulator", "target", "weights", "resc_weights"]
        output = {0: pd.DataFrame([["A", "C", 1, 1], ["B", "D", 1, 1], ["A", "D", 1, 1]], columns=cols),
                  1: pd.DataFrame([["B", "D", 1, 1], ["C", "D", 1, 1], ["A", "F", 1, 1]], columns=cols)}

        targets = pd.Index(["A", "B", "C", "D", "E", "F"], name="target")
        regulators = pd.Index(["A", "B", "C"], name="regulator")

        correct_output_0 = pd.DataFrame(0., index=targets, columns=regulators)
        correct_output_0.loc["C", "A"] = 1.
        correct_output_0.loc["D", "B"] = 1.
        correct_output_0.loc["D", "A"] = 1.

        out_0 = amusr_regression._format_weights(output[0], "weights", targets, regulators)
        pdt.assert_frame_equal(out_0, correct_output_0)

    def test_process_output_all(self):
        cols = ["regulator", "target", "weights", "resc_weights"]
        output = {0: pd.DataFrame([["A", "C", 1, 1], ["B", "D", 1, 1], ["A", "D", 1, 1]], columns=cols),
                  1: pd.DataFrame([["B", "D", 1, 1], ["C", "D", 1, 1], ["A", "F", 1, 1]], columns=cols)}

        targets = pd.Index(["A", "B", "C", "D", "E", "F"], name="target")
        regulators = pd.Index(["A", "B", "C"], name="regulator")

        correct_output_0 = pd.DataFrame(0., index=targets, columns=regulators)
        correct_output_0.loc["C", "A"] = 1.
        correct_output_0.loc["D", "B"] = 1.
        correct_output_0.loc["D", "A"] = 1.

        correct_output_1 = pd.DataFrame(0., index=targets, columns=regulators)
        correct_output_1.loc["D", "B"] = 1.
        correct_output_1.loc["D", "C"] = 1.
        correct_output_1.loc["F", "A"] = 1.

        runner = amusr_regression.AMuSR_regression([pd.DataFrame()], [pd.DataFrame()], genes=targets.tolist(),
                                                   tfs=regulators.tolist())
        runner.n_tasks = 2

        weights, resc_weights = runner.pileup_data([output])

        pdt.assert_frame_equal(weights[0], correct_output_0)
        pdt.assert_frame_equal(resc_weights[0], correct_output_0)
        pdt.assert_frame_equal(weights[1], correct_output_1)
        pdt.assert_frame_equal(resc_weights[1], correct_output_1)
