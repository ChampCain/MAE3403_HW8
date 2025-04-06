# region imports
import math
from Calc_state import *
from UnitConversions import UnitConverter as UC
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy as dc


# these imports are necessary for drawing a matplot lib graph on my GUI
# no simple widget for this exists in QT Designer, so I have to add the widget in code.
# endregion

# region class definitions
class rankineModel():
    def __init__(self):
        '''
        Constructor for rankine power cycle data (in the Model-View-Controller design pattern).
        This class is for storing data only. The Controller class should update the model depending
        on input from the user. The View class should display the model data depending on the desired output.
        :param p_low: the low pressure isobar for the cycle in kPa
        :param p_high: the high pressure isobar for the cycle in kPa
        :param t_high: optional temperature for State1 (turbine inlet) in degrees C
        :param eff_turbine: isentropic efficiency of the turbine
        :param name: a convenient name
        '''
        self.p_low = None  # in kPa
        self.p_high = None  # in kPa
        self.t_high = None  # in C
        self.name = None
        self.efficiency = None
        self.turbine_eff = None
        self.turbine_work = None
        self.pump_work = None
        self.heat_added = None
        self.steam = Steam_SI()  # Instantiate a steam object for calculating state
        # Initialize the states as stateProps objects
        self.state1 = stateProps()
        self.state2s = stateProps()
        self.state2 = stateProps()
        self.state3 = stateProps()
        self.state4 = stateProps()
        self.SI = True  # If False, convert pressures from psi to kPa and T from F to C for inputs
        # Data for plotting
        self.satLiqPlotData = StateDataForPlotting()
        self.satVapPlotData = StateDataForPlotting()
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()


class rankineView():
    def __init__(self):
        """
        Empty constructor by design
        """
        pass

    def setWidgets(self, *args):
        # Create class variables for the input widgets
        self.rb_SI, self.le_PHigh, self.le_PLow, self.le_TurbineInletCondition, \
            self.rdo_Quality, self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis, \
            self.chk_logX, self.chk_logY = args[0]

        # Create class variables for the display widgets
        self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow, self.lbl_SatPropHigh, \
            self.lbl_TurbineInletCondition, self.lbl_H1, self.lbl_H1Units, self.lbl_H2, \
            self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units, self.lbl_H4, self.lbl_H4Units, \
            self.lbl_TurbineWork, self.lbl_TurbineWorkUnits, self.lbl_PumpWork, \
            self.lbl_PumpWorkUnits, self.lbl_HeatAdded, self.lbl_HeatAddedUnits, \
            self.lbl_ThermalEfficiency, self.canvas, self.figure, self.ax = args[1]

    def selectQualityOrTHigh(self, Model=None):
        """
        Action to take when selecting one of the radio buttons for Quality or THigh
        """
        if Model is None:
            return

        SI = self.rb_SI.isChecked()
        if self.rdo_Quality.isChecked():
            self.le_TurbineInletCondition.setText("1.0")
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            try:
                # Get PHigh from the line edit and convert to kPa
                PHigh_text = self.le_PHigh.text()
                if not PHigh_text:
                    return

                PCF = 100 if SI else UC.psi_to_kpa  # Convert to kPa
                PHigh = float(PHigh_text) * PCF

                # Get saturation temperature at PHigh
                satPHigh = Model.steam.getsatProps_p(PHigh)
                if satPHigh is None:
                    print(f"Could not get saturation properties for P={PHigh} kPa")
                    return

                Tsat = satPHigh.tsat

                # Convert temperature if needed and update the line edit
                if not SI:
                    Tsat = UC.C_to_F(Tsat)
                self.le_TurbineInletCondition.setText("{:0.2f}".format(Tsat))
                self.le_TurbineInletCondition.setEnabled(True)
            except ValueError:
                print("Invalid pressure value entered")
            except Exception as e:
                print(f"Error in selectQualityOrTHigh: {e}")

        x = self.rdo_Quality.isChecked()
        self.lbl_TurbineInletCondition.setText(
            ("Turbine Inlet: {}{} =".format('x' if x else 'THigh', '' if x else ('(C)' if SI else '(F)'))))

    def setNewPHigh(self, Model=None):
        """
        Update the saturated properties display when P High changes
        """
        if Model is None:
            return

        try:
            SI = self.rb_SI.isChecked()
            # Convert input to kPa (steam tables use kPa)
            PCF = 100 if SI else UC.psi_to_kpa  # bar->kPa or psi->kPa
            PHigh_text = self.le_PHigh.text()

            if not PHigh_text:
                return

            PHigh = float(PHigh_text) * PCF

            # Validate pressure range (typical steam tables go up to ~220 bar)
            if PHigh > 22000:  # 220 bar = 22000 kPa
                print(f"Warning: Pressure {PHigh} kPa is too high for steam tables")
                return

            # Get saturation properties
            satPHigh = Model.steam.getsatProps_p(PHigh)
            if satPHigh is None:
                print(f"Could not get saturation properties for P={PHigh} kPa")
                return

            # Update saturation properties display
            self.lbl_SatPropHigh.setText(satPHigh.getTextOutput(SI=SI))

            # If in T High mode, update the turbine inlet temperature
            if self.rdo_THigh.isChecked():
                Tsat = satPHigh.tsat
                if not SI:
                    Tsat = UC.C_to_F(Tsat)
                self.le_TurbineInletCondition.setText(f"{Tsat:.2f}")
        except ValueError:
            print("Invalid pressure value entered")
        except Exception as e:
            print(f"Error in setNewPHigh: {e}")

    def setNewPLow(self, Model=None):
        """
        Update the saturated properties display when P Low changes
        """
        if Model is None:
            return

        try:
            SI = self.rb_SI.isChecked()
            # Convert input to kPa (steam tables use kPa)
            PCF = 100 if SI else UC.psi_to_kpa  # bar->kPa or psi->kPa
            PLow_text = self.le_PLow.text()

            if not PLow_text:
                return

            PLow = float(PLow_text) * PCF

            # Get saturation properties
            satPLow = Model.steam.getsatProps_p(PLow)
            if satPLow is None:
                print(f"Could not get saturation properties for P={PLow} kPa")
                return

            # Update saturation properties display
            self.lbl_SatPropLow.setText(satPLow.getTextOutput(SI=SI))
        except ValueError:
            print("Invalid pressure value entered")
        except Exception as e:
            print(f"Error in setNewPLow: {e}")

    def outputToGUI(self, Model=None):
        if Model is None or Model.state1 is None:  # means the cycle has not been evaluated yet
            return

        # Update the line edits and labels
        HCF = 1 if Model.SI else UC.kJperkg_to_BTUperlb  # Enthalpy conversion factor
        self.lbl_H1.setText("{:0.2f}".format(Model.state1.h * HCF))
        self.lbl_H2.setText("{:0.2f}".format(Model.state2.h * HCF))
        self.lbl_H3.setText("{:0.2f}".format(Model.state3.h * HCF))
        self.lbl_H4.setText("{:0.2f}".format(Model.state4.h * HCF))
        self.lbl_TurbineWork.setText("{:0.2f}".format(Model.turbine_work * HCF))
        self.lbl_PumpWork.setText("{:0.2f}".format(Model.pump_work * HCF))
        self.lbl_HeatAdded.setText("{:0.2f}".format(Model.heat_added * HCF))
        self.lbl_ThermalEfficiency.setText("{:0.2f}".format(Model.efficiency))

        # Update saturation properties
        satPropsLow = Model.steam.getsatProps_p(p=Model.p_low)
        satPropsHigh = Model.steam.getsatProps_p(p=Model.p_high)
        self.lbl_SatPropLow.setText(satPropsLow.getTextOutput(SI=Model.SI))
        self.lbl_SatPropHigh.setText(satPropsHigh.getTextOutput(SI=Model.SI))

        # Update the plot
        self.plot_cycle_XY(Model=Model)

    def updateUnits(self, Model=None):
        """
        Updates the units on the GUI to match choice of SI or English
        """
        if Model is None:
            return

        # Update the outputs first
        self.outputToGUI(Model=Model)

        # Update pressure labels
        SI = Model.SI
        self.lbl_PHigh.setText("P High (bar)" if SI else "P High (psi)")
        self.lbl_PLow.setText("P Low (bar)" if SI else "P Low (psi)")

        # Update THigh if in THigh mode
        if not self.rdo_Quality.isChecked():
            self.selectQualityOrTHigh(Model)

        # Update energy units
        energyUnits = "kJ/kg" if SI else "BTU/lb"
        self.lbl_H1Units.setText(energyUnits)
        self.lbl_H2Units.setText(energyUnits)
        self.lbl_H3Units.setText(energyUnits)
        self.lbl_H4Units.setText(energyUnits)
        self.lbl_TurbineWorkUnits.setText(energyUnits)
        self.lbl_PumpWorkUnits.setText(energyUnits)
        self.lbl_HeatAddedUnits.setText(energyUnits)

        # Update the plot with new units
        self.plot_cycle_XY(Model=Model)

    def print_summary(self, Model=None):
        """
        Prints to CLI.
        """
        if Model is None:
            return

        if Model.efficiency is None:
            Model.calc_efficiency()

        print('Cycle Summary for: ', Model.name)
        print('\tEfficiency: {:0.3f}%'.format(Model.efficiency))
        print('\tTurbine Eff:  {:0.2f}'.format(Model.turbine_eff))
        print('\tTurbine Work: {:0.3f} kJ/kg'.format(Model.turbine_work))
        print('\tPump Work: {:0.3f} kJ/kg'.format(Model.pump_work))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(Model.heat_added))
        Model.state1.print()
        Model.state2.print()
        Model.state3.print()
        Model.state4.print()

    def plot_cycle_TS(self, axObj=None, Model=None):
        """
        Plots the Rankine cycle on T-S coordinates along with the vapor dome
        """
        if Model is None:
            return

        SI = Model.SI
        steam = Model.steam

        # Step 1&2: Build data for saturated liquid and vapor lines
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt('sat_water_table.txt',
                                                          skiprows=1, unpack=True)
        ax = plt.subplot() if axObj is None else axObj

        # Unit conversion factors
        hCF = 1 if SI else UC.kJperkg_to_BTUperlb
        pCF = 1 if SI else UC.kpa_to_psi
        sCF = 1 if SI else UC.kJperkgK_to_BTUperlbR
        vCF = 1 if SI else UC.kgperm3_to_lbperft3

        # Apply unit conversions
        sfs *= sCF
        sgs *= sCF
        hfs *= hCF
        hgs *= hCF
        vfs *= vCF
        vgs *= vCF
        ps *= pCF
        ts = [t if SI else UC.C_to_F(t) for t in ts]

        # Plot saturated liquid and vapor lines
        ax.plot(sfs, ts, color='blue')
        ax.plot(sgs, ts, color='red')

        # Step 3: Line between state3 and saturated liquid at p_high
        st3p = steam.getState(Model.p_high, x=0)  # saturated liquid state at p_high
        svals = np.linspace(Model.state3.s, st3p.s, 20)
        tvals = np.linspace(Model.state3.T, st3p.T, 20)
        line3 = np.column_stack([svals, tvals])

        # Step 4: Line between saturated liquid and vapor at p_high
        sat_pHigh = steam.getState(Model.p_high, x=1.0)
        st1 = Model.state1
        svals2p = np.linspace(st3p.s, sat_pHigh.s, 20)
        tvals2p = [st3p.T for _ in range(20)]
        line4 = np.column_stack([svals2p, tvals2p])

        if st1.T > sat_pHigh.T:  # Add data points for superheated region
            svals_sh = np.linspace(sat_pHigh.s, st1.s, 20)
            tvals_sh = np.array([steam.getState(Model.p_high, s=ss).T for ss in svals_sh])
            line4 = np.append(line4, np.column_stack([svals_sh, tvals_sh]), axis=0)

        # Step 5: Line between state1 and state2
        svals = np.linspace(Model.state1.s, Model.state2.s, 20)
        tvals = np.linspace(Model.state1.T, Model.state2.T, 20)
        line5 = np.column_stack([svals, tvals])

        # Step 6: Line between state2 and state3
        svals = np.linspace(Model.state2.s, Model.state3.s, 20)
        tvals = np.array([Model.state2.T for _ in range(20)])
        line6 = np.column_stack([svals, tvals])

        # Step 7: Combine lines for upper curve
        topLine = np.append(line3, line4, axis=0)
        topLine = np.append(topLine, line5, axis=0)
        xvals = topLine[:, 0]
        y1 = topLine[:, 1]
        y2 = [Model.state3.T for _ in xvals]

        # Convert units if English
        if not SI:
            xvals *= UC.kJperkgK_to_BTUperlbR
            y1 = [UC.C_to_F(t) for t in y1]
            y2 = [UC.C_to_F(t) for t in y2]

        # Plot the cycle lines
        ax.plot(xvals, y1, color='darkgreen')
        ax.plot(xvals, y2, color='black')

        # Plot state points
        if SI:
            ax.plot(Model.state1.s, Model.state1.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s, Model.state2.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s, Model.state3.T, marker='o', markeredgecolor='k', markerfacecolor='w')
        else:
            ax.plot(Model.state1.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state1.T),
                    marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state2.T),
                    marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state3.T),
                    marker='o', markeredgecolor='k', markerfacecolor='w')

        # Set labels and title
        tempUnits = r'$\left(^oC\right)$' if SI else r'$\left(^oF\right)$'
        entropyUnits = r'$\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r'$\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        ax.set_xlabel(r's ' + entropyUnits, fontsize=18)
        ax.set_ylabel(r'T ' + tempUnits, fontsize=18)
        ax.set_title(Model.name)
        ax.grid(visible='both', alpha=0.5)
        ax.tick_params(axis='both', direction='in', labelsize=18)

        # Set axis limits
        sMin = min(sfs)
        sMax = max(sgs)
        ax.set_xlim(sMin, sMax)

        tMin = min(ts)
        tMax = max(max(ts), st1.T)
        ax.set_ylim(tMin, tMax * 1.05)

        if axObj is None:  # Show plot if not being displayed on a figure
            plt.show()

    def plot_cycle_XY(self, Model=None):
        """
        Plots any two thermodynamic properties on X and Y axes
        """
        if Model is None:
            return

        ax = self.ax
        X = self.cmb_XAxis.currentText()
        Y = self.cmb_YAxis.currentText()
        logx = self.chk_logX.isChecked()
        logy = self.chk_logY.isChecked()
        SI = Model.SI

        if X == Y:
            return

        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if ax is None:
            ax = plt.subplot()
            QTPlotting = False

        ax.clear()
        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')

        # Get data for plotting
        YF = Model.satLiqPlotData.getDataCol(Y, SI=SI)
        YG = Model.satVapPlotData.getDataCol(Y, SI=SI)
        XF = Model.satLiqPlotData.getDataCol(X, SI=SI)
        XG = Model.satVapPlotData.getDataCol(X, SI=SI)

        # Plot the vapor dome
        ax.plot(XF, YF, color='b')
        ax.plot(XG, YG, color='r')

        # Plot the upper and lower curves
        ax.plot(Model.lowerCurve.getDataCol(X, SI=SI),
                Model.lowerCurve.getDataCol(Y, SI=SI), color='k')
        ax.plot(Model.upperCurve.getDataCol(X, SI=SI),
                Model.upperCurve.getDataCol(Y, SI=SI), color='g')

        # Add axis labels
        ax.set_ylabel(Model.lowerCurve.getAxisLabel(Y, SI=SI),
                      fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel(Model.lowerCurve.getAxisLabel(X, SI=SI),
                      fontsize='large' if QTPlotting else 'medium')

        # Add title
        Model.name = 'Rankine Cycle - ' + Model.state1.region + ' at Turbine Inlet'
        ax.set_title(Model.name, fontsize='large' if QTPlotting else 'medium')

        # Format tick marks
        ax.tick_params(axis='both', which='both', direction='in',
                       top=True, right=True,
                       labelsize='large' if QTPlotting else 'medium')

        # Plot state points
        ax.plot(Model.state1.getVal(X, SI=SI), Model.state1.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state2.getVal(X, SI=SI), Model.state2.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state3.getVal(X, SI=SI), Model.state3.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(Model.state4.getVal(X, SI=SI), Model.state4.getVal(Y, SI=SI),
                marker='o', markerfacecolor='w', markeredgecolor='k')

        # Set axis limits
        xmin = min(min(XF), min(XG), min(Model.upperCurve.getDataCol(X, SI=SI)),
                   min(Model.lowerCurve.getDataCol(X, SI=SI)))
        xmax = max(max(XF), max(XG), max(Model.upperCurve.getDataCol(X, SI=SI)),
                   max(Model.lowerCurve.getDataCol(X, SI=SI)))
        ymin = min(min(YF), min(YG), min(Model.upperCurve.getDataCol(Y, SI=SI)),
                   min(Model.lowerCurve.getDataCol(Y, SI=SI)))
        ymax = max(max(YF), max(YG), max(Model.upperCurve.getDataCol(Y, SI=SI)),
                   max(Model.lowerCurve.getDataCol(Y, SI=SI))) * 1.1

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        if not QTPlotting:
            plt.show()
        else:
            self.canvas.draw()


class rankineController():
    def __init__(self, *args):
        """
        Create rankineModel object. The rankineController class updates the model
        based on user input and updates the rankineView as well
        """
        self.Model = rankineModel()
        self.View = rankineView()
        self.IW = args[0]  # input widgets
        self.DW = args[1]  # display widgets
        self.View.setWidgets(self.IW, self.DW)
        self.buildVaporDomeData()

    def updateModel(self):
        """
        Read input widgets and update the model
        """
        self.Model.SI = self.View.rb_SI.isChecked()

        # Convert input to proper units (kPa for internal calculations)
        PCF = 100 if self.Model.SI else UC.psi_to_kpa  # bar->kPa or psi->kPa

        try:
            self.Model.p_high = float(self.View.le_PHigh.text()) * PCF
            self.Model.p_low = float(self.View.le_PLow.text()) * PCF
            T = float(self.View.le_TurbineInletCondition.text())
            self.Model.t_high = None if self.View.rdo_Quality.isChecked() else (
                T if self.Model.SI else UC.F_to_C(T))
            self.Model.turbine_eff = float(self.View.le_TurbineEff.text())

            # Perform calculations
            self.calc_efficiency()
            self.updateView()
        except Exception as e:
            print(f"Error in updateModel: {e}")

    def updateUnits(self):
        """
        Update units throughout the GUI
        """
        self.Model.SI = self.View.rb_SI.isChecked()
        self.View.updateUnits(Model=self.Model)

    def selectQualityOrTHigh(self):
        """
        Handle quality vs T High selection
        """
        self.View.selectQualityOrTHigh(self.Model)

    def setNewPHigh(self):
        """
        Handle P High changes
        """
        self.View.setNewPHigh(self.Model)

    def setNewPLow(self):
        """
        Handle P Low changes
        """
        self.View.setNewPLow(self.Model)

    def calc_efficiency(self):
        """
        Calculate the Rankine cycle efficiency
        """
        steam = self.Model.steam

        # State 1: turbine inlet (p_high, t_high) superheated or saturated vapor
        if self.Model.t_high is None:
            self.Model.state1 = steam.getState(P=self.Model.p_high, x=1.0, name='Turbine Inlet')
        else:
            self.Model.state1 = steam.getState(P=self.Model.p_high, T=self.Model.t_high, name='Turbine Inlet')

        # State 2s: turbine exit (p_low, s=s_turbine inlet) isentropic
        self.Model.state2s = steam.getState(P=self.Model.p_low, s=self.Model.state1.s, name="Turbine Exit")

        # State 2: actual turbine exit considering efficiency
        if self.Model.turbine_eff < 1.0:
            h2 = self.Model.state1.h - self.Model.turbine_eff * (self.Model.state1.h - self.Model.state2s.h)
            self.Model.state2 = steam.getState(P=self.Model.p_low, h=h2, name="Turbine Exit")
        else:
            self.Model.state2 = self.Model.state2s

        # State 3: pump inlet (p_low, x=0) saturated liquid
        self.Model.state3 = steam.getState(P=self.Model.p_low, x=0, name='Pump Inlet')

        # State 4: pump exit (p_high, s=s_pump_inlet)
        self.Model.state4 = steam.getState(P=self.Model.p_high, s=self.Model.state3.s, name='Pump Exit')

        # Calculate work and heat terms
        self.Model.turbine_work = self.Model.state1.h - self.Model.state2.h
        self.Model.pump_work = self.Model.state4.h - self.Model.state3.h
        self.Model.heat_added = self.Model.state1.h - self.Model.state4.h
        self.Model.efficiency = 100.0 * (self.Model.turbine_work - self.Model.pump_work) / self.Model.heat_added

        return self.Model.efficiency

    def updateView(self):
        """
        Update the view with current model data
        """
        self.buildDataForPlotting()
        self.View.outputToGUI(Model=self.Model)

    def setRankine(self, p_low=8, p_high=8000, t_high=None, eff_turbine=1.0, name='Rankine Cycle'):
        """
        Set model values for rankine power cycle
        """
        self.Model.p_low = p_low  # in kPa
        self.Model.p_high = p_high  # in kPa
        self.Model.t_high = t_high  # in C
        self.Model.name = name
        self.Model.efficiency = None
        self.Model.turbine_eff = eff_turbine
        self.Model.turbine_work = 0
        self.Model.pump_work = 0
        self.Model.heat_added = 0
        self.Model.state1 = None  # entrance to turbine
        self.Model.state2s = None  # entrance to condenser (isentropic turbine)
        self.Model.state2 = None  # entrance to condenser (non-isentropic turbine)
        self.Model.state3 = None  # entrance to pump (saturated liquid at plow)
        self.Model.state4 = None  # entrance to boiler (isentropic)

    def print_summary(self):
        """
        Print cycle summary to console
        """
        self.View.print_summary(Model=self.Model)

    def buildVaporDomeData(self, nPoints=500):
        """
        Build the vapor dome from just above the triple point up to the critical point
        """
        steam = self.Model.steam
        tp = triplePt_PT()
        cp = criticalPt_PT()

        steam.state.p = cp.p
        steam.state.t = cp.t
        steam.calcState_1Phase()
        critProps = dc(steam.state)

        P = np.logspace(math.log10(tp.p * 1.001), math.log10(cp.p * 0.99), nPoints)

        for p in P:
            sat = steam.getsatProps_p(p)
            self.Model.satLiqPlotData.addPt((sat.tsat, p, sat.uf, sat.hf, sat.sf, sat.vf))
            self.Model.satVapPlotData.addPt((sat.tsat, p, sat.uf, sat.hg, sat.sg, sat.vg))

        self.Model.satLiqPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))
        self.Model.satVapPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))

    def buildDataForPlotting(self):
        """
        Create data for plotting the Rankine cycle
        """
        # Clear old data
        self.Model.upperCurve.clear()
        self.Model.lowerCurve.clear()

        steam = self.Model.steam
        satPLow = steam.getsatProps_p(self.Model.p_low)
        satPHigh = steam.getsatProps_p(self.Model.p_high)

        # Build upper curve (processes 3-4-1-2)
        # Region 3-4: Pump (isentropic compression)
        nPts = 15
        DeltaP = (satPHigh.psat - satPLow.psat)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=(satPLow.psat + z * DeltaP), s=satPLow.sf)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Region 4-1: Heat addition in boiler
        # First from T4 to T5 where T5 is the saturated liquid at p_high
        T4 = state.t
        T5 = satPHigh.tsat
        DeltaT = (T5 - T4)
        nPts = 20
        P = satPHigh.psat

        for n in range(nPts - 1):
            z = n * 1.0 / (nPts - 2)
            T = T4 + z * DeltaT
            if T < T5:
                state = steam.getState(P=P, T=T)
                self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Then through the two-phase region
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(satPHigh.psat, x=z)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Finally through superheated region if needed
        if self.Model.state1.t > (satPHigh.tsat + 1):
            T6 = satPHigh.tsat
            DeltaT = self.Model.state1.t - T6
            for n in range(nPts):
                z = n * 1.0 / (nPts - 1)
                if z > 0:
                    state = steam.getState(satPHigh.psat, T=T6 + z * DeltaT)
                    self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Region 1-2: Turbine expansion
        s1 = self.Model.state1.s
        s2 = self.Model.state2.s
        P1 = self.Model.state1.p
        P2 = self.Model.state2.p
        Deltas = s2 - s1
        DeltaP = P2 - P1

        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=P1 + z * DeltaP, s=s1 + z * Deltas)
            self.Model.upperCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Build lower curve (process 2-3: Condensation)
        x2 = self.Model.state2.x
        state = self.Model.state2

        # Account for superheated vapor at turbine exit
        if state.t > satPLow.tsat:
            nPts = 20
            DeltaT = (state.t - satPLow.tsat) / nPts
            self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

            for n in range(nPts):
                t = self.Model.state2.t - n * DeltaT
                if t > satPLow.tsat:
                    state = steam.getState(P=satPLow.psat, T=t)
                    self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

        # Through the two-phase region
        nPts = len(self.Model.upperCurve.t)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=satPLow.psat, x=(1.0 - z) * x2)
            self.Model.lowerCurve.addPt((state.t, state.p, state.u, state.h, state.s, state.v))

    def updatePlot(self):
        """
        Update the plot with current data
        """
        self.View.plot_cycle_XY(Model=self.Model)


# endregion

# region function definitions
def main():
    RC = rankineController()
    RC.setRankine(8, 8000, t_high=500, eff_turbine=0.9,
                  name='Rankine Cycle - Superheated at turbine inlet')
    eff = RC.calc_efficiency()
    print(eff)
    RC.print_summary()
    RC.View.plot_cycle_TS(Model=RC.Model)


# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion