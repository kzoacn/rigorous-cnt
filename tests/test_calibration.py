"""Checks for explicit spectral calibration; this is separate from core sampling."""

import unittest
from dataclasses import replace
from sage.all import QQ,ZZ,RealBallField,PolynomialRing,NumberField
from rigorous_cnt import (NumberFieldContext,RayConditions,calibrate_walk,
                          algorithm2,RandomBits,ResourceLimit)
from rigorous_cnt.numerics import endpoints


R=PolynomialRing(QQ,'t'); t=R.gen()


class CalibrationTests(unittest.TestCase):
    def test_finite_class_group_exact_mixing_comparison(self):
        K=NumberField(t*t+5,'a'); ctx=NumberFieldContext(K); ray=RayConditions(ctx)
        walk=calibrate_walk(ctx,ray,max_prime_bound=76)
        certificate=walk.certificate
        self.assertEqual(certificate.cyclic_orders,(2,))
        self.assertEqual(certificate.low_characters,1)
        # For B=19, the two-state walk has nontrivial eigenvalue -2/3.
        exact_l1=(QQ(2)/3)**walk.steps
        self.assertLessEqual(exact_l1,walk.mixing_l1_bound)
        self.assertLessEqual(walk.mixing_l1_bound,QQ(1)/200)
        result=algorithm2(ctx,K.ideal(1),walk,ray=ray,rng=RandomBits(11))
        self.assertEqual(result.report['mixing_guarantee'],'verified_spectral_bound')

    def test_real_quadratic_frequency_coverage_and_bound(self):
        K=NumberField(t*t-2,'a'); ctx=NumberFieldContext(K); ray=RayConditions(ctx)
        walk=calibrate_walk(ctx,ray,max_prime_bound=76)
        c=walk.certificate
        RB=RealBallField(512)
        unit_length=RB(2).sqrt()*(1+RB(2).sqrt()).log()
        lo,hi=endpoints(unit_length*c.cutoff_radius)
        self.assertEqual(lo.floor(),hi.floor())
        self.assertEqual(c.low_characters,2*lo.floor())
        s=QQ(1)/4; sprime=QQ(1)/11200
        independently_computed=(2*RB(c.volume_upper)/sprime*(
            RB(c.spectral_upper)**(2*c.steps)+(-2*(RB(c.cutoff_radius)*s)**2).exp())).sqrt()
        self.assertLessEqual(endpoints(independently_computed)[1],QQ(1)/200)
        result=algorithm2(ctx,K.ideal(3),walk,ray=ray,rng=RandomBits(3))
        self.assertIn(result.element,K.ideal(3))

    def test_nontrivial_ray_with_partial_real_signs(self):
        K=NumberField(t*t-2,'a'); ctx=NumberFieldContext(K)
        ray=RayConditions(ctx,K.ideal(3),real_places=[0],tau=-1)
        walk=calibrate_walk(ctx,ray,max_prime_bound=304,max_frequency_vectors=10000)
        result=algorithm2(ctx,K.ideal(2),walk,ray=ray,rng=RandomBits(7))
        self.assertTrue(ray.contains(result.element))
        self.assertEqual(result.report['mixing_guarantee'],'verified_spectral_bound')

    def test_mixed_cubic_spectral_calibration(self):
        K=NumberField(t*t*t-2,'a'); ctx=NumberFieldContext(K); ray=RayConditions(ctx)
        walk=calibrate_walk(ctx,ray,max_prime_bound=152)
        self.assertGreater(walk.certificate.low_characters,0)
        self.assertLessEqual(walk.mixing_l1_bound,QQ(1)/200)

    def test_totally_complex_quartic_with_nonzero_unit_rank(self):
        K=NumberField(t**4+1,'a'); ctx=NumberFieldContext(K); ray=RayConditions(ctx)
        walk=calibrate_walk(ctx,ray,max_prime_bound=304)
        result=algorithm2(ctx,K.ideal(3),walk,ray=ray,rng=RandomBits(21))
        self.assertEqual(ctx.weights,(2,2))
        self.assertEqual(sum(result.log_distortion),0)
        self.assertEqual(result.distortion[0]**2*result.distortion[1]**2,1)
        self.assertIn(result.element,K.ideal(3))

    def test_certificate_is_bound_to_parameters_and_ray(self):
        K=NumberField(t*t+5,'a'); ctx=NumberFieldContext(K); ray=RayConditions(ctx)
        walk=calibrate_walk(ctx,ray,max_prime_bound=76)
        with self.assertRaises(ValueError):
            algorithm2(ctx,K.ideal(1),replace(walk,steps=walk.steps+1),ray=ray)
        with self.assertRaises(ValueError):
            algorithm2(ctx,K.ideal(1),walk,ray=RayConditions(ctx,K.ideal(3)))
        with self.assertRaises(ValueError):
            algorithm2(ctx,K.ideal(1),walk,ray=ray,subgroup_oracle=lambda _:True)

    def test_calibration_resource_limit(self):
        K=NumberField(t*t-2,'a'); ctx=NumberFieldContext(K)
        with self.assertRaises(ResourceLimit):
            calibrate_walk(ctx,max_frequency_vectors=1)


if __name__=='__main__': unittest.main()
