// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::LiquidPhaseThreeC
 */
#ifndef DUMUX_LIQUID_THREEC_PHASE_HH
#define DUMUX_LIQUID_THREEC_PHASE_HH

#include <cassert>
#include <limits>

#include <dune/common/exceptions.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/binarycoefficients/h2o_constant.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {
/*!
 * \ingroup Fluidsystems
 * \brief Policy for the H2O-N2 fluid system
 */
template<bool fastButSimplifiedRelations = true>
struct liquidp3cDefaultPolicy
{
    static constexpr  bool useH2ODensityAsLiquidMixtureDensity() { return fastButSimplifiedRelations; }
};



/*!
 * \ingroup Fluidsystems
 * \brief A liquid phase consisting of a two components,
 *        a main component and a conservative tracer component
 */
template <class Scalar, class MainComponent, class SecondComponent, class ThirdComponent, class Policy = liquidp3cDefaultPolicy<>>
class LiquidPhaseThreeC
: public Base<Scalar, LiquidPhaseThreeC<Scalar, MainComponent, SecondComponent, ThirdComponent> >
{
    using ThisType = LiquidPhaseThreeC<Scalar, MainComponent, SecondComponent, ThirdComponent,Policy>;
    //using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    //file only for h2o + one componant. better re-write own liquidDiffCoeff() function
	//using BinaryCoefficients = BinaryCoeff::H2O_Component<Scalar, SecondComponent>;//??
    using BinaryCoefficientsSecondComponent = BinaryCoeff::H2O_Component<Scalar, SecondComponent>;
    using BinaryCoefficientsThirdComponent  = BinaryCoeff::H2O_Component<Scalar, ThirdComponent>;

public:
    using ParameterCache = NullParameterCache;

    static constexpr int numPhases = 1; //!< Number of phases in the fluid system
    static constexpr int numComponents = 3; //!< Number of components in the fluid system

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the only phase

    static constexpr int comp0Idx = 0; //!< index of the frist component
    static constexpr int comp1Idx = 1; //!< index of the second component
    static constexpr int comp2Idx = 2; //!< index of the 3rd component
    static constexpr int mainCompIdx = comp0Idx; //!< index of the main component
    static constexpr int secondCompIdx = comp1Idx; //!< index of the secondary component
    static constexpr int thirdCompIdx = comp2Idx; //!< index of the secondary component

    /*!
    * \brief Initialize the fluid system's static parameters generically
    */
    static void init() {}

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return IOName::liquidPhase(); }

    /*!
     * \brief Returns whether the fluids are miscible
     * \note There is only one phase, so miscibility makes no sense
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    { 
		std::string componentName_;
		switch(compIdx)
		{
			case comp0Idx:
			{
				componentName_ = MainComponent::name();
				break;
			}
			case comp1Idx:
			{
				componentName_ = SecondComponent::name();
				break;
			}
			case comp2Idx:
			{
				componentName_ = ThirdComponent::name();
				break;
			}
			default:
			{
				DUNE_THROW(Dune::InvalidStateException, "liquidphase3c::componentName: compIdx "<<compIdx<<" not recognised");
			}				
		}		
		return componentName_; 
	}

    /*!
     * \brief A human readable name for the fluid system.
     */
    static std::string name()
    { return "LiquidPhaseThreeC"; }

    /*!
     * \brief Returns whether the fluid is gaseous
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if only a single component is involved. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx = 0)
    { return true; }

    /*!
     * \brief Returns true if the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return MainComponent::liquidIsCompressible(); }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static bool isIdealGas(int phaseIdx = 0)
    { return false; /* we're a liquid! */ }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx)
    {  
		Scalar molarMass_;
		switch(compIdx)
		{
			case comp0Idx:
			{
				molarMass_ = MainComponent::molarMass();
				break;
			}
			default:
			{
				DUNE_THROW(Dune::InvalidStateException, "liquidphase3c::molarMass: compIdx "<<compIdx<<" not recognised");
			}
		}		
		return molarMass_; 
	}

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of the main component
     */
    static Scalar criticalTemperature()
    {  return MainComponent::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of the main component
     */
    static Scalar criticalPressure()
    {  return MainComponent::criticalPressure(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at the main component's triple point.
     */
    static Scalar tripleTemperature()
    {  return MainComponent::tripleTemperature(); }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at the main component's triple point.
     */
    static Scalar triplePressure()
    { return MainComponent::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the main component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    {  return MainComponent::vaporPressure(T); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the phase at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return MainComponent::liquidDensity(temperature, pressure); }

    using Base<Scalar, ThisType>::density;
    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the phase at a given pressure and temperature.
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx = 0)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // See: Eq. (7) in Class et al. (2002a)
        // This assumes each gas molecule displaces exactly one
        // molecule in the liquid.
        const Scalar pureComponentMolarDensity = MainComponent::liquidMolarDensity(T, p);
		//return pureComponentMolarDensity *fluidState.moleFraction(phase0Idx, mainCompIdx);
		// assume pure water
		
            if (Policy::useH2ODensityAsLiquidMixtureDensity())
			{
                // assume pure water
                return MainComponent::liquidDensity(T, p);
			}
            else
            {
        return pureComponentMolarDensity
               * (MainComponent::molarMass()*fluidState.moleFraction(phase0Idx, mainCompIdx)
                  + SecondComponent::molarMass()*fluidState.moleFraction(phase0Idx, secondCompIdx)
                  + ThirdComponent::molarMass()*fluidState.moleFraction(phase0Idx, thirdCompIdx));
			}
				  
				  
    }

    using Base<Scalar, ThisType>::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // assume pure component or that each gas molecule displaces exactly one
        // molecule in the liquid.
        return MainComponent::liquidMolarDensity(T, p);
    }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return MainComponent::liquidPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return MainComponent::liquidEnthalpy(temperature, pressure); }

    using Base<Scalar, ThisType>::enthalpy;
    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const int phaseIdx)
    {
        return enthalpy(fluidState.temperature(phaseIdx),
                        fluidState.pressure(phaseIdx));
    }

    /*!
    * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in the specified phase
    * \param fluidState The fluid state
    * \param phaseIdx The index of the phase
    * \param componentIdx The index of the component
    */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int componentIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (componentIdx == mainCompIdx)
            return MainComponent::liquidEnthalpy(T, p);
        else if (componentIdx == secondCompIdx)
            return SecondComponent::liquidEnthalpy(T, p);
        else if (componentIdx == thirdCompIdx)
            return ThirdComponent::liquidEnthalpy(T, p);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
    }


    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return MainComponent::liquidInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return MainComponent::liquidViscosity(temperature, pressure); }

    using Base<Scalar, ThisType>::viscosity;
    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        return viscosity(fluidState.temperature(phaseIdx),
                         fluidState.pressure(phaseIdx));
    }

    using Base<Scalar, ThisType>::fugacityCoefficient;
    /*!
     * \copybrief Base::fugacityCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    using Base<Scalar, ThisType>::diffusionCoefficient;
    /*!
     * \copybrief Base::diffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
	 * @see the equation from h2o_constant
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException, "Not applicable: Diffusion coefficients");
    }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    /*!
     * \copybrief Base::binaryDiffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the background (i think?) component to consider
     * \param compJIdx The index of the component to consider
	 *
	 * Binary diffusion coefficient => diffusion of compJIdx in compIIdx I believe
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState, int phaseIdx, int compIIdx, int compJIdx)
    {
        assert(phaseIdx == compIIdx);
		
        if (compJIdx == secondCompIdx)
            return BinaryCoefficientsSecondComponent::liquidDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        else if (compJIdx == thirdCompIdx)
            return BinaryCoefficientsThirdComponent::liquidDiffCoeff(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        return -1;
    }
	

    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return MainComponent::liquidThermalConductivity(temperature, pressure); }

    using Base<Scalar, ThisType>::thermalConductivity;
    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        return thermalConductivity(fluidState.temperature(phaseIdx),
                                   fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return MainComponent::liquidHeatCapacity(temperature, pressure); }

    using Base<Scalar, ThisType>::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        return heatCapacity(fluidState.temperature(phaseIdx),
                            fluidState.pressure(phaseIdx));
    }
	
};

} // namespace FluidSystems

} // namespace Dumux

#endif
