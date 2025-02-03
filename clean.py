import os


if __name__ == "__main__":
   #get pathway to current directory
   pathway = os.getcwd()
   os.listdir(pathway)

   #files to delete with extensinos
   files_delete = ['.bai', 'output.tsv']

   #go through files in directory
   for file in os.listdir(pathway):
      for extension in files_delete:
         #delete files we want to delete
         if extension in file:
            curr_dir = os.getcwd()
            delete_pathway = os.path.join(curr_dir, file)
            os.remove(delete_pathway)
