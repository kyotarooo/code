class Book:
    def __init__(self,title,author):#title: タイトル　author: 著者
        self.title = title
        self.author = author

    def show_info(self):
        return f"{self.title}: {self.author}"#↑ここまで(1)
    
book1 = Book("走れメロス","太宰治")
book2 = Book("こころ","夏目漱石")#↑ここまで(2)

print(f'{book1.show_info()}')
print(f'{book2.show_info()}')#↑ここまで(3)


