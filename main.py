from tkinter import Tk
from tkinter import Label
from tkinter import ttk



if __name__  == '__main__':
 
    root = Tk()     # создаем корневой объект - окно
    root.title("Курсовая работа ПиАЭ")     # устанавливаем заголовок окна
    root.geometry("1280x960")    # устанавливаем размеры окна
    
    label = Label(text="Hello METANIT.COM") # создаем текстовую метку
    label.pack()    # размещаем метку в окне
    
    tabControl = ttk.Notebook(root)
    tab1 = ttk.Frame(tabControl)
    tab2 = ttk.Frame(tabControl)
    tab3 = ttk.Frame(tabControl)
    tab4 = ttk.Frame(tabControl)
    tabControl.add(tab1, text ='Синтез')
    tabControl.add(tab2, text ='Графики')
    tabControl.add(tab3, text ='Таблицы')
    tabControl.add(tab4, text ='О приложении')
    tabControl.pack(expand = 1, fill ="both")
    
    
    root.mainloop()
    
   
