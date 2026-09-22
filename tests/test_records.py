"""Independent checks of saved end-to-end outputs, without rerunning sampling."""

import json
import unittest
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,NumberField
from rigorous_cnt import NumberFieldContext,FactorBase,CompactElement,verify_sunit_record
from rigorous_cnt.compact import combine_compact


RESULTS=Path(__file__).resolve().parents[1]/'examples'/'results'


class RecordedResultsTests(unittest.TestCase):
    def test_saved_results_are_recertified(self):
        for name in ('real_units.json','class2_sunits.json'):
            data=json.loads((RESULTS/name).read_text())
            data['report']={'complete_generation':'this flag is deliberately untrusted'}
            result=verify_sunit_record(data)
            self.assertTrue(result['verified'])
            self.assertFalse(result['stored_flags_trusted'])

    def test_full_rank_proper_subgroup_record_is_rejected(self):
        data=json.loads((RESULTS/'real_units.json').read_text())
        data['exponent_matrix'][0]=[str(2*ZZ(x)) for x in data['exponent_matrix'][0]]
        # This is still a valid full-rank unit subgroup, but has index two.
        result=verify_sunit_record(data)
        self.assertFalse(result['verified'])
        self.assertEqual(result['target_index_interval'],(2,2))

    def test_materialized_results_against_sage(self):
        for name in ('real_units.json','class2_sunits.json'):
            data=json.loads((RESULTS/name).read_text())
            R=PolynomialRing(QQ,'t')
            K=NumberField(R([QQ(x) for x in data['defining_polynomial']]),'a')
            ctx=NumberFieldContext(K)
            def ideal(rows):return K.ideal([K([QQ(x) for x in row]) for row in rows])
            working=FactorBase(ctx,[ideal(rows) for rows in data['working_prime_ideals']])
            requested=[ideal(rows) for rows in data['requested_prime_ideals']]
            source=[CompactElement(ctx,[(K([QQ(x) for x in term['coefficients']]),ZZ(term['exponent']))
                                        for term in factors]) for factors in data['source_elements']]
            generators=[combine_compact(ctx,source,[ZZ(x) for x in row]).materialize(working)
                        for row in data['exponent_matrix']]
            reference=K.S_unit_group(S=requested,proof=True)
            self.assertEqual(len(generators),1)
            self.assertEqual(abs(ZZ(reference.log(generators[0])[1])),1)


if __name__=='__main__':unittest.main()
