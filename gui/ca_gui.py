from tkinter import *
from PIL import ImageTk, Image
import os

import print_parameter_files

root = Tk()
root.title("Ca simulation GUI")
home =os.path.expanduser("~")
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
ca_image_file = os.path.join(__location__, 'images/milch.png')
loop_image_file = os.path.join(__location__, 'images/loop.png')
ca_img = ImageTk.PhotoImage(Image.open(ca_image_file))
loop_img = ImageTk.PhotoImage(Image.open(loop_image_file).resize((20, 20), Image.ANTIALIAS))

# Set parameter main frame
parameter_frame = LabelFrame(root, text="Parameters", padx=5, pady=5)
parameter_frame.pack(padx=10, pady=10)

# Set Frame for IP3R parameters
ip3r_parameter_frame = LabelFrame(parameter_frame, text="IP3R", padx=5, pady=5)
ip3r_parameter_frame.pack(padx=50, pady=10)

# Set Frame for global Ca dynamics
ca_parameter_frame = LabelFrame(parameter_frame, text="Ca2+", padx=5, pady=5)
ca_parameter_frame.pack(padx=50, pady=10)

# Set Frame for global Inhibitory dynamics
inhib_parameter_frame = LabelFrame(parameter_frame, text="Inhibition", padx=5, pady=5)
inhib_parameter_frame.pack(padx=50, pady=10)


# Set Entries and respective districtions (labels)
ip3r_label1 =  Label(ip3r_parameter_frame, text="Num. Sites")
ip3r_label2 =  Label(ip3r_parameter_frame, text="Num. Channels")
ip3r_label3 =  Label(ip3r_parameter_frame, text="tau")
ip3r_label4 =  Label(ip3r_parameter_frame, text="mu")
ip3r_label5 =  Label(ip3r_parameter_frame, text="D")

ip3r_entry1 = Entry(ip3r_parameter_frame, text="Num. Sites")
ip3r_entry2 = Entry(ip3r_parameter_frame, text="Num. Channels")
ip3r_entry3 = Entry(ip3r_parameter_frame, text="tau")
ip3r_entry4 = Entry(ip3r_parameter_frame, text="mu")
ip3r_entry5 = Entry(ip3r_parameter_frame, text="D")

def add_loop_settings(parent):
    ent1 = Entry(parent, text="Start")
    ent2 = Entry(parent, text="Stop")
    ent3 = Entry(parent, text="Steps")
    ent1.grid(row=0, column=1)
    ent2.grid(row=0, column=2)
    ent3.grid(row=0, column=3)

ip3r_button1 = Button(ip3r_parameter_frame, image=loop_img, height=20, width=20, command=lambda: add_loop_settings(ip3r_parameter_frame))
ip3r_button2 = Button(ip3r_parameter_frame, image=loop_img, height=20, width=20, command=add_loop_settings)
ip3r_button3 = Button(ip3r_parameter_frame, image=loop_img, height=20, width=20, command=add_loop_settings)
ip3r_button4 = Button(ip3r_parameter_frame, image=loop_img, height=20, width=20, command=add_loop_settings)
ip3r_button5 = Button(ip3r_parameter_frame, image=loop_img, height=20, width=20, command=add_loop_settings)

ca_label1 = Label(ca_parameter_frame, text="tau_ca")
ca_entry1 = Entry(ca_parameter_frame, text="tau_ca")
ca_label2 = Label(ca_parameter_frame, text="alpha")
ca_entry2 = Entry(ca_parameter_frame, text="alpha")

inhib_label1 = Label(inhib_parameter_frame, text="tau_inhib")
inhib_entry1 = Entry(inhib_parameter_frame, text="tau_inhib")
inhib_label2 = Label(inhib_parameter_frame, text="beta")
inhib_entry2 = Entry(inhib_parameter_frame, text="beta")

ip3r_entry1.grid(row=0, column=1)
ip3r_entry2.grid(row=1, column=1)
ip3r_entry3.grid(row=2, column=1)
ip3r_entry4.grid(row=3, column=1)
ip3r_entry5.grid(row=4, column=1)
ip3r_label1.grid(row=0, column=0)
ip3r_label2.grid(row=1, column=0)
ip3r_label3.grid(row=2, column=0)
ip3r_label4.grid(row=3, column=0)
ip3r_label5.grid(row=4, column=0)
ip3r_button1.grid(row=0, column=2)
ip3r_button2.grid(row=1, column=2)
ip3r_button3.grid(row=2, column=2)
ip3r_button4.grid(row=3, column=2)
ip3r_button5.grid(row=4, column=2)

ca_entry1.grid(row=0, column=1)
ca_entry2.grid(row=1, column=1)
ca_label1.grid(row=0, column=0)
ca_label2.grid(row=1, column=0)

inhib_entry1.grid(row=0, column=1)
inhib_entry2.grid(row=1, column=1)
inhib_label1.grid(row=0, column=0)
inhib_label2.grid(row=1, column=0)

def save_parameters():
    N = ip3r_entry1.get()
    M = ip3r_entry2.get()
    tau_x = ip3r_entry3.get()
    mu_x = ip3r_entry4.get()
    D_x = ip3r_entry5.get()
    tau_ca = ca_entry1.get()
    alpha_ca = ca_entry2.get()
    tau_inhib = inhib_entry1.get()
    beta_inhib = inhib_entry2.get()
    print_parameter_files.print_parameter_files(N, M, tau_x, mu_x, D_x, tau_ca, alpha_ca, tau_inhib, beta_inhib)

save_parameter_button = Button(parameter_frame, text="Save parameters", command=save_parameters)
save_parameter_button.pack()


root.mainloop()