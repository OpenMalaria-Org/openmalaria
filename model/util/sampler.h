/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef UTIL_SAMPLER
#define UTIL_SAMPLER

#include <limits>
#include <optional>

#include "Global.h"
#include "util/random.h"
#include "schema/util.h"

namespace OM { namespace util {
    
        /** A normal sample, which can be turned into various log-normal samples.
     * 
     * Allows generation of correlated log-normal samples with different sigma.
     */
    class NormalSample {
    public:
        NormalSample() : x(numeric_limits<double>::signaling_NaN()) {}
        
        /// convert sample to N(mu,sigma)
        double asNormal( double mu, double sigma ) const;

        /// convert sample to lnN(mu,sigma)
        double asLognormal( double mu, double sigma ) const;
        
        static NormalSample generate(LocalRng& rng);
        
        /// Generate a sample correlated with base.
        /// 
        /// @param base sample correlated against
        /// @param correlation expected correlation between these parameters
        ///     (note: if this is used to sample a log-normal, then the
        ///     correlation is on the log-scale)
        /// @param factor should equal sqrt(1 - correlation^2); may be cached
        static NormalSample generate_correlated(NormalSample base, double correlation, double factor, LocalRng& rng);
        
    private:
        NormalSample( double variate ) : x(variate) {}
        
        double x;   // variate (sampled from N(0,1))
    };
    
    /** Sampler for normal values */
    class NormalSampler {
    public:
        NormalSampler() :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters such that samples taken are:
         * X ~ N( m, s² )
         * 
         * @param m Mean of sampled variates
         * @param s Square-root of variance of sampled variates
         */
        void setParams( double m, double s );

        /** As above, using an XML element. */
        void setParams( const scnXml::SampledValueN& elt );
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;

        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample sample) const{
            return sample.asNormal( mu, sigma );
        }
        
        /// Return mu / mean of distribution
        inline double getMu() const{
            return mu;
        }

        /// Return sigma / standard deviation of distribution
        inline double getSigma() const{
            return sigma;
        }
        
    private:
        double mu, sigma;
    };
    
    class Sampler
    {
    public:
        virtual ~Sampler() = default;
        virtual std::unique_ptr<Sampler> clone() const = 0;
        virtual double sample(LocalRng& rng) const = 0;
        virtual double sample(NormalSample sample) const{ throw std::runtime_error("sample(NormalSample) not implemented for this distribution");}
        virtual double mean() const = 0;
        virtual double cdf(double x) const { throw std::runtime_error("cdf() not implemented for this distribution"); }
    };
    
    /** Sampler for log-normal values */
    class LognormalSampler : public Sampler {
    public:
        LognormalSampler() = default;
        LognormalSampler& operator=(const LognormalSampler&) = default;

        std::unique_ptr<Sampler> clone() const override {
            return std::make_unique<LognormalSampler>(*this);
        }

        /**
         * Configure the log-normal distribution from a target mean and coefficient of variation.
         *
         * For CV > 0, this computes the underlying normal parameters μ and σ such that
         *  E[X] = exp(μ + σ²/2) = mean
         *  CV   = √(exp(σ²) − 1)
         *
         * If CV == 0, the distribution collapses to a point mass:
         *  - σ is set to 0
         *  - μ is set to ln(mean), or –∞ if mean == 0
         *  - sampling and CDF methods treat X == exp(μ)
         *
         * @param mean  Desired mean of X (must be ≥ 0; if mean == 0 then CV must be 0).
         * @param CV    Desired coefficient of variation (must be ≥ 0).
         * @throws util::xml_scenario_error  If mean < 0 or CV < 0, or if mean > 0 but CV == 0 is invalid.
         */
        static unique_ptr<util::LognormalSampler> fromMeanCV( double mean, double CV, std::optional<double> truncate = std::nullopt);

        /**
         * Configure the log-normal distribution from a target mean and variance.
         *
         * For variance > 0, this computes the equivalent CV = variance/mean and then
         * derives μ and σ via the same formulas as in setMeanCV.
         *
         * If variance == 0, the distribution collapses to a point mass:
         *  - σ is set to 0
         *  - μ is set to ln(mean)
         *  - sampling and CDF methods treat X == exp(μ)
         *
         * @param mean      Desired mean of X (must be > 0).
         * @param variance  Desired variance of X (must be ≥ 0).
         * @throws util::xml_scenario_error  If mean ≤ 0 or variance < 0.
         */
        static unique_ptr<util::LognormalSampler> fromMeanVariance( double mean, double variance, std::optional<double> truncate = std::nullopt);

        /** Scale the mean (i.e. multiply by a scalar).
         * 
         * Note that sigma (when specified via CV) is independent of the mean,
         * so one can safely multiply the mean without affecting CV. */
        void scaleMean( double scalar );
        
        /** Get the mean. */
        double mean() const override;

        /** Sample a value. */
        double sample(LocalRng& rng) const override;
        
        /** Create a log-normal sample from an existing normal sample. */
        double sample(NormalSample sample) const override {
            return sample.asLognormal( mu, sigma );
        }

        /**
         * Compute the cumulative distribution function (CDF) of this log-normal sampler.
         *
         * Special cases:
         *  - If x ≤ 0, returns 0.
         *  - If sigma == 0, the distribution collapses to a point mass,
         * so returns 1 if log(x) ≥ mu, else 0.
         *  - Otherwise uses gsl_cdf_lognormal_P(x, mu, sigma).
         *
         * @param x  The point at which to evaluate the CDF (must be a real number).
         * @return   P(X ≤ x), where X ~ LogNormal(mu, sigma²).
         */
        double cdf(double x) const override;
        
    private:
        // log-space parameters
        double mu       = std::numeric_limits<double>::signaling_NaN(); 
        double sigma    = std::numeric_limits<double>::signaling_NaN();
        double CV       = std::numeric_limits<double>::signaling_NaN();
        std::optional<double> truncate = std::nullopt;
    };

    /** Sampler for gamma values */
    class GammaSampler : public Sampler {
    public:
        GammaSampler() = default;
        GammaSampler& operator=(const GammaSampler&) = default;

        std::unique_ptr<Sampler> clone() const override {
            return std::make_unique<GammaSampler>(*this);
        }

        /**
         * Configure the gamma distribution from a target mean and coefficient of variation.
         *
         * Computes the shape (k) and scale (θ) parameters so that
         *  mean = k * θ
         *  CV   = 1 / √k
         *
         * If CV == 0, the distribution collapses to a point mass; in this case
         *  - k and θ are set to NaN
         *  - sampling and CDF methods treat this as a degenerate distribution.
         *
         * @param mean  Desired mean of the distribution (must be > 0).
         * @param CV    Desired coefficient of variation (must be ≥ 0).
         * @throws util::xml_scenario_error  If mean ≤ 0 or CV < 0.
         */
        static unique_ptr<util::GammaSampler> fromMeanCV( double mean, double CV, std::optional<double> truncate = std::nullopt);

        /**
         * Configure the gamma distribution from a target mean and variance.
         *
         * Computes the shape (k) and scale (θ) parameters so that
         *  mean     = k * θ
         *  variance = k * θ²
         *
         * If variance == 0, the distribution collapses to a point mass; in this case
         *  - k and θ are set to NaN
         *  - sampling and CDF methods treat this as a degenerate distribution at `mean`.
         *
         * @param mean      Desired mean of the distribution (must be > 0).
         * @param variance  Desired variance of the distribution (must be ≥ 0).
         * @throws util::xml_scenario_error  If mean ≤ 0 or variance < 0.
         */
        static unique_ptr<util::GammaSampler> fromMeanVariance( double mean, double variance, std::optional<double> truncate = std::nullopt);
        
        /** Get the mean. */
        double mean() const override;

        /** Sample a value. */
        double sample(LocalRng& rng) const override;
        
        /**
         * Compute the cumulative distribution function (CDF) of this gamma sampler.
         *
         * Special cases:
         *  - If x ≤ 0, returns 0.
         *  - If CV == 0 or variance == 0, the distribution collapses to a point mass,
         * so returns 1 if x ≥ mu, else 0.
         *  - Otherwise uses gsl_cdf_gamma_P(x, k, theta) where
         *      k     = shape parameter
         *      theta = scale parameter
         *
         * @param x  The point at which to evaluate the CDF (must be a real number).
         * @return   P(X ≤ x), where X ~ Gamma(k, θ).
         */
        double cdf(double x) const override;
        
    private:
        double mu       = std::numeric_limits<double>::signaling_NaN();
        double k        = std::numeric_limits<double>::signaling_NaN();
        double theta    = std::numeric_limits<double>::signaling_NaN();
        double variance = std::numeric_limits<double>::signaling_NaN();
        double CV       = std::numeric_limits<double>::signaling_NaN();
        std::optional<double> truncate = std::nullopt;
    };

    /** Create a new <SamplerT> from a xml snippet. */
    template <typename SamplerT>
    inline unique_ptr<SamplerT> createSampler(double mean, const scnXml::SampledValueCV& elt)
    {
        // 1.  Compile‑time identify which concrete sampler we are handling ─────
        constexpr bool isLognormal = std::is_same_v<SamplerT, LognormalSampler>;
        constexpr bool isGamma     = std::is_same_v<SamplerT, GammaSampler>;
        constexpr const char* expectedDistr =
            isLognormal ? "lognormal" :
            isGamma     ? "gamma"
                        : nullptr;

        const std::string& d = elt.getDistr();

        // 2. Special handling for the “const” alias (lognormal with CV = 0)
        if( d == "const" )
        {
            if constexpr (!isLognormal) // only allowed for lognormal
                throw util::xml_scenario_error("\"distr=\"const\": expected \"lognormal\"");

            if( elt.getCV().present() && elt.getCV().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"CV\" must be zero or omitted when distr=\"const\" or is omitted" );
            if( elt.getVariance().present() && elt.getVariance().get() != 0.0 )
                throw util::xml_scenario_error( "\"distr="+d+"\": attribute \"variance\" must be zero or omitted when distr=\"const\" or is omitted" );
            
            return SamplerT::fromMeanCV(mean, 0.0);
        }

        // 3. Check that XML 'distr' matches the compile‑time sampler type
        if (expectedDistr == nullptr)
            throw util::xml_scenario_error("createSampler<SamplerT>: expected SamplerT to be a valid Sampler");

        if (d != expectedDistr)
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"" + expectedDistr + "\" or \"const\"");

        // 4. Exactly one of CV / variance must be supplied
        const bool hasCV  = elt.getCV().present();
        const bool hasVar = elt.getVariance().present();

        if (hasCV == hasVar)   // both present or both absent
            throw util::xml_scenario_error("\"distr=" + d + "\": exactly one of attributes \"CV\" or \"variance\" must be specified");

        // 5. Check for the optional truncation parameter
        std::optional<double> truncate = elt.getTruncate().present()
            ? std::optional<double>{elt.getTruncate().get()}
            : std::nullopt;

        // 6. Build the sampler
        if(hasCV)
            return SamplerT::fromMeanCV(mean, elt.getCV().get(), truncate);
        else
            return SamplerT::fromMeanVariance(mean, elt.getVariance().get(), truncate);
    }

    /** Create a Sampler from a SampledValueCV xml snippet. */
    inline std::unique_ptr<util::Sampler> createSampler(double mean, const scnXml::SampledValueCV& elt)
    {
        const std::string& d = elt.getDistr();

        if (d == "const" || d == "lognormal")
            return createSampler<util::LognormalSampler>(mean, elt);
        else if (d == "gamma")
            return createSampler<util::GammaSampler>(mean,  elt);
        else
            throw util::xml_scenario_error("\"distr=" + d + "\": expected \"const\", \"lognormal\" or \"gamma\"");
    }

    /** Create a LognormalSampler from a SampledValueLN xml snippet. */
    inline std::unique_ptr<LognormalSampler> createSampler(const scnXml::SampledValueLN& elt)
    {
        return createSampler<LognormalSampler>(elt.getMean(), static_cast<const scnXml::SampledValueCV&>(elt));
    }

    /** Sampler for the Beta distribution.
     *
     * Input may be alpha and beta or mean and variance. As a special case,
     * variance zero is also supported, meaning return the mean. */
    class BetaSampler {
    public:
        BetaSampler() :
            a( numeric_limits<double>::signaling_NaN() ),
            b( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters: alpha, beta.
         * 
         * As a special case, when beta=0, alpha is interpreted as the mean and
         * is returned without sampling. */
        inline void setParams( double alpha, double beta ){
            a=alpha;
            b=beta;
        }
        
        /** Set parameters: alpha and beta are calculated such that the mean
         * and variance are as given. Variance may be zero, in which case the
         * mean is returned directly without sampling. */
        void setParamsMV( double mean, double variance );
        
        /** Set parameters from an XML element. */
        void setParams( const scnXml::BetaMeanSample& elt ){
            setParamsMV( elt.getMean(), elt.getVariance() );
        }
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;
        
    private:
        //Note: if b is 0, then alpha is mean. Otherwise a and b are
        //the expected (α,β) parameters to the beta distribution.
        double a, b;
    };
    
    /** Sampler for the Weibull distribution. */
    class WeibullSampler {
    public:
        WeibullSampler() :
            scale( numeric_limits<double>::signaling_NaN() ),
            shape( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters from scale, shape. */
        inline void setScaleShape( double lambda, double k ){
            scale = lambda;
            shape = k;
        }
        
        /** Set parameters from an XML element. */
        inline void setParams( const scnXml::WeibullSample& elt ){
            setScaleShape( elt.getScale(), elt.getShape() );
        }
        
        /** Sample a value. */
        inline double sample(LocalRng& rng) const{
            return rng.weibull( scale, shape );
        }
        
    private:
        double scale, shape;    // λ, k
    };
    
} }
#endif
