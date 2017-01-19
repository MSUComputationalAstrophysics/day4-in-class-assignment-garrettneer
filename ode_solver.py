"""
Title: Day 4 PCA: ODE Solver
Author: Garrett Neer
Date: 1-18-17
Description: Solves ODE (spring mass system) with various methods. 
Euler, predictor-corrector, runge-kutta... 

"""
from math import pi,sin,cos,fabs
from ROOT import TGraph, TMultiGraph,TF1,TCanvas,TLegend,kBlue,kGreen,kBlack,kRed,gStyle

def calculate_energy(position,velocity):
	return 0.5*position**2 + 0.5*velocity**2

def euler(position_list,velocity_list, time_step):
	current_position = position_list[len(position_list)-1]
	current_velocity = velocity_list[len(velocity_list)-1]
	position_list.append(current_position+current_velocity*time_step)
	velocity_list.append(current_velocity-current_position*time_step)

def picard(position_list,velocity_list,time_step):
	current_position = position_list[len(position_list)-2]
	current_velocity = velocity_list[len(velocity_list)-2]
	position_list[len(position_list)-1] = current_position + (current_velocity + velocity_list[len(velocity_list)-1])*time_step/2
	velocity_list[len(velocity_list)-1] = current_velocity - (current_position + position_list[len(position_list)-1])*time_step/2

def predictor_corrector(position_list,velocity_list,time_step):
	euler(position_list,velocity_list,time_step)
	picard(position_list,velocity_list,time_step)

def runge_kutta(position_list,velocity_list,time_step):
	current_position = position_list[len(position_list)-1]
	current_velocity = velocity_list[len(velocity_list)-1]
	x1 = current_position
	v1 = current_velocity
	a1 = -current_position
	x2 = current_position + time_step*v1/2
	v2 = current_velocity + time_step*a1/2
	a2 = -x2
	x3 = current_position + time_step*v2/2
	v3 = current_velocity + time_step*a2/2
	a3 = -x3
	x4 = current_position + time_step*v3
	v4 = current_velocity + time_step*a3
	a4 = -x4
	new_position = current_position + time_step*(v1+2*v2+2*v3+v4)/6
	new_velocity = current_velocity + time_step*(a1+2*a2+2*a3+a4)/6
	position_list.append(new_position)
	velocity_list.append(new_velocity)

def midpoint(position_list,velocity_list,time_step):
	current_position = position_list[-1]
	current_velocity = velocity_list[-1]
	current_acceleration = -current_position
	midpoint_position = current_position+current_velocity*time_step/2
	midpoint_velocity = current_velocity+time_step*current_acceleration/2
	midpoint_acceleration = -midpoint_position
	position_list.append(current_position+midpoint_velocity*time_step)
	velocity_list.append(current_velocity+midpoint_acceleration*time_step)

time_step_list = [0.1*pi,0.01*pi,0.001*pi]
omega = 1
constant = 1
initial_time = 0
final_time = 4*pi
initial_energy_list = []
final_energy_list = []

from array import array
euler_graph_list = []
predictor_corrector_graph_list = []
runge_kutta_graph_list = []
midpoint_graph_list = []
marker_list = [24,25,26]
color_list = [kBlue, kGreen+1,kBlack]
def solve_ode(ode_method, graph_list):
	for time_step in time_step_list:
		current_time = initial_time
		x_list = [0.0]
		position_list = [0.0]
		velocity_list = [1.0]
		initial_energy_list.append(calculate_energy(position_list[len(position_list)-1],velocity_list[len(velocity_list)-1]))
		while current_time <= final_time:
			ode_method(position_list,velocity_list,time_step)
			current_time += time_step
			x_list.append(current_time)

		final_energy_list.append(calculate_energy(position_list[len(position_list)-1],velocity_list[len(velocity_list)-1]))
		x_array = array("d",x_list)
		position_array = array("d",position_list)
		graph_list.append(TGraph(len(x_list),x_array,position_array))



solve_ode(euler,euler_graph_list)
solve_ode(predictor_corrector,predictor_corrector_graph_list)
solve_ode(runge_kutta,runge_kutta_graph_list)
solve_ode(midpoint,midpoint_graph_list)


energy_change_list = []
for i in range(0,len(initial_energy_list)):
	energy_change_list.append(fabs(final_energy_list[i]-initial_energy_list[i])/initial_energy_list[i])


energy_change_marker_list = [24,25,26,27]
energy_change_color_list = [kRed,kBlack,kBlue,kGreen+1]
energy_change_graph_list = []
for i in range(0,4):
	energy_change_graph_list.append(TGraph(len(time_step_list),array("d",time_step_list),array("d",energy_change_list[i*3:i*3+3])))
	energy_change_graph_list[i].SetMarkerStyle(energy_change_marker_list[i])
	energy_change_graph_list[i].SetMarkerColor(energy_change_color_list[i])
	energy_change_graph_list[i].SetLineColor(energy_change_color_list[i])
	energy_change_graph_list[i].SetLineWidth(2)

energy_change_mg = TMultiGraph()
for i in range(0,4):
	energy_change_mg.Add(energy_change_graph_list[i])

euler_multigraph = TMultiGraph()
predictor_corrector_multigraph = TMultiGraph()
runge_kutta_multigraph = TMultiGraph()
midpoint_multigraph = TMultiGraph()
for i in range(0,len(time_step_list)):
	euler_graph_list[i].SetMarkerStyle(marker_list[i])
	euler_graph_list[i].SetMarkerColor(color_list[i])
	midpoint_graph_list[i].SetMarkerStyle(marker_list[i])
	midpoint_graph_list[i].SetMarkerColor(color_list[i])
	runge_kutta_graph_list[i].SetMarkerStyle(marker_list[i])
	runge_kutta_graph_list[i].SetMarkerColor(color_list[i])
	predictor_corrector_graph_list[i].SetMarkerStyle(marker_list[i])
	predictor_corrector_graph_list[i].SetMarkerColor(color_list[i])
	euler_multigraph.Add(euler_graph_list[i])
	midpoint_multigraph.Add(midpoint_graph_list[i])
	runge_kutta_multigraph.Add(runge_kutta_graph_list[i])
	predictor_corrector_multigraph.Add(predictor_corrector_graph_list[i])

analytic_position = TF1("","sin(x)",initial_time,final_time)
analytic_position.SetTitle("")
analytic_position.GetXaxis().SetTitle("time (s)")
analytic_position.GetXaxis().SetTitle("position (m)")

legend_label_list = ["time step: 0.1#pi","time step: 0.01#pi","time step: 0.001#pi"]
gStyle.SetLegendBorderSize(0)
c1 = TCanvas()
euler_multigraph.Draw("AP")
euler_multigraph.SetTitle("Euler Method")
euler_multigraph.GetXaxis().SetTitle("time (s)")
euler_multigraph.GetYaxis().SetTitle("position (m)")
legend = TLegend(0.69,0.69,0.89,0.89)
for i in range(0,len(euler_graph_list)):
	legend.AddEntry(euler_graph_list[i],legend_label_list[i],"p")
legend.Draw()
analytic_position.Draw("sames")
c1.SaveAs("euler.png")

c2 = TCanvas()
predictor_corrector_multigraph.Draw("AP")
predictor_corrector_multigraph.SetTitle("Predictor-Corrector Method")
predictor_corrector_multigraph.GetXaxis().SetTitle("time (s)")
predictor_corrector_multigraph.GetYaxis().SetTitle("position (m)")
legend = TLegend(0.69,0.69,0.89,0.89)
for i in range(0,len(predictor_corrector_graph_list)):
	legend.AddEntry(predictor_corrector_graph_list[i],legend_label_list[i],"p")
legend.Draw()
analytic_position.Draw("sames")
c2.SaveAs("predictor_corrector.png")

energy_change_label_list = ["Euler","Predictor-Corrector","Runge-Kutta","Midpoint"]
c3 = TCanvas()
c3.SetLogx()
c3.SetLogy()
energy_change_mg.Draw("APL")
energy_change_mg.GetXaxis().SetTitle("time step (s)")
energy_change_mg.GetYaxis().SetTitle("#epsilon")
legend = TLegend(0.54,0.21,0.84,0.40)
for i in range(0,len(energy_change_graph_list)):
	legend.AddEntry(energy_change_graph_list[i],energy_change_label_list[i],"p")
legend.Draw()

c3.SaveAs("energy_change.png")

c4 = TCanvas()
runge_kutta_multigraph.Draw("AP")
runge_kutta_multigraph.SetTitle("4th order Runge-Kutta Method")
runge_kutta_multigraph.GetXaxis().SetTitle("time (s)")
runge_kutta_multigraph.GetYaxis().SetTitle("position (m)")
legend = TLegend(0.69,0.69,0.89,0.89)
for i in range(0,len(runge_kutta_graph_list)):
	legend.AddEntry(runge_kutta_graph_list[i],legend_label_list[i],"p")
legend.Draw()
analytic_position.Draw("sames")
c1.SaveAs("runge_kutta.png")

c5 = TCanvas()
midpoint_multigraph.Draw("AP")
midpoint_multigraph.SetTitle("Midpoint Method")
midpoint_multigraph.GetXaxis().SetTitle("time (s)")
midpoint_multigraph.GetYaxis().SetTitle("position (m)")
legend = TLegend(0.69,0.69,0.89,0.89)
for i in range(0,len(midpoint_graph_list)):
	legend.AddEntry(midpoint_graph_list[i],legend_label_list[i],"p")
legend.Draw()
analytic_position.Draw("sames")
c5.SaveAs("midpoint.png")

raw_input()
