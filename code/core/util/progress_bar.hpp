#ifndef _progress_bar_hpp_
#define _progress_bar_hpp_

class progress_bar {
  public:
    progress_bar(int _size, real _v0, real _v1) : size(_size), v0(_v0), v1(_v1), perc(0) {

      std::cout << "[0%";
      for(int i = 0; i < size / 2 - 2; ++i)
        std::cout << " ";
      std::cout << "|50%";
      for(int i = 0; i < (size - size / 2) - 8; ++i)
        std::cout << " ";
      std::cout << "100%]\n[";
    }
    void update(real value) {
      real new_perc = (value - v0) / (v1 - v0);
      int dx = std::min(static_cast<int>(new_perc * size), size) - perc;

      for(int i = 0; i < dx; ++i)
        std::cout << ".";
      perc += dx;
    }

    void cancel() {
      if(perc < size) {
        std::cout << 'X';
        ++perc;
        for(int i = perc; i < size; ++i)
          std::cout << ' ';
      }
      perc = size;
    }

    void done() {
      int dx = size - perc;

      for(int i = 0; i < dx; ++i)
        std::cout << ".";
      perc = size;
      std::cout << "]\n";
    }

  private:
    real v0, v1;
    int size;
    int perc;
};

#endif // _progress_bar_hpp_
