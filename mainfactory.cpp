#include <iostream>
#include <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>
#include <chrono>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){

    auto &factory = Factory::Instance();

    // Read data from main arg
    std::ifstream file;
    if(argc < 2){
        std::cerr << "Error: no filename given for data as arg" << std::endl;
        return 1;
    }
    file.open(argv[1]);
    if(!file.is_open()){
        std::cerr << "Error: " << argv[1] << " file does not exist" <<
            std::endl;
        return 1;
    }
    std::string str, str2;
    std::getline(file, str);
    std::istringstream iss(str);

    std::vector<double> v;
  
    while(std::getline(iss, str2, ',')){
        double val = ::atof(str2.c_str());
        v.push_back(val);
    }

    file.close();
    Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        v.data(), v.size()); // TODO: meglio con conservative resize?

    double mu0, lambda, alpha0, beta0;
    //std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
    //std::cin >> mu0 >> lambda >> alpha0 >> beta0; // 5.0 0.1 2.0 2.0
    // HypersFixedNNIG hy(mu0, lambda, alpha0, beta0);

    double totalmass;
    //std::cout << "Insert total mass value:" << std::endl; 
    //std::cin >> totalmass; //1.0
    // DirichletMixture mix(totalmass);

    int n_aux = 3;
    //std::cout << "Insert number of auxiliary blocks:" << std::endl;
    //std::cin >> n_aux;

    //std::ofstream file;
    //file.open("csv/data.csv");
    //for(auto &d : data){
    //    file << d << ",";
    //}
    //file << std::endl;
    //file.close();
 
    HypersFixedNNIG hy(5.0, 1.0, 2.0, 2.0); // mu0, lambda, alpha0, beta0
    DirichletMixture mix(1); // total mass

    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
        data, n_aux, mix, hy);
  
    // Run sampler
    BaseCollector *f;
    if(argc < 3){
        std::cerr << "Error: need file collector type " <<
            "(\"file\" or \"memory\") as arg" << std::endl;
        return 1;
    }

    std::string collector(argv[2]);
    if(collector == "file"){
        std::string filename;
        if(argc < 4){
            // Use default name
            filename = "collector.bin";
        }
        else {
            std::string filename = argv[2];
            if(argc > 4){
                std::cout << "Warning: unused extra args present" << std::endl;
            }
        }
        f = new FileCollector(filename);
    }

    else if(collector == "memory"){
        if(argc > 3){
            std::cout << "Warning: unused extra args present" << std::endl;
        }
        f = new MemoryCollector();
    }

    else {
        std::cerr << "Error: collector type arg must be \"file\" or \"memory\""
            << std::endl;
        return 1;
    }
  
    sampler.run(f);

    // Density and clustering stuff stuff
    double temp = 0.0;
    double step = 0.05;
    double upp_bnd = 10.0;
    std::vector<double> v_temp;
    while(temp <= upp_bnd){
        v_temp.push_back(temp);
        temp += step;
    }
    Eigen::VectorXd grid = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        v_temp.data(), v_temp.size()); 

    sampler.eval_density(grid, f);
    sampler.write_density_to_file("csv/density_ex.csv");
    unsigned int i_cap = sampler.cluster_estimate(f);
    std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    sampler.write_final_clustering_to_file();
    sampler.write_best_clustering_to_file();

    return 0;
}















//! A simple example of use of the generic object Factory template class
#include "Factory.hpp"
#include <iostream>
// an abstract class
class Abstract
{
public:
  virtual void whoAmI()=0;
  virtual ~Abstract()=default;
};

//! Some conclrete classes
/*!@{ */
class Derived1: public Abstract
{
public:
  void whoAmI()
  {
    std::cout<<" I am Derived1"<<std::endl;
  }
};

class Derived2: public Abstract
{
public:
  void whoAmI()
  {
    std::cout<<" I am Derived2"<<std::endl;
  };
};

class Derived3: public Abstract
{
public:
  void whoAmI()
  {
    std::cout<<" I am Derived3"<<std::endl;
  };
};
/*!@} */

using Identifier=int;
using Builder=std::function<std::unique_ptr<Abstract>()>;
using MyObjectFactory=GenericFactory::Factory<Abstract,Identifier,Builder>;

// Defining the builders for the different concrete objects
Builder build1=[]{return std::make_unique<Derived1>();};
Builder build2=[]{return std::make_unique<Derived2>();};
Builder build3=[]{return std::make_unique<Derived3>();};


  
//! Normally this is done elsewhere, but this is only a test
/*! 
  Filling up the factory with the builders 
 */
void loadFactory()
{
  auto & factory=MyObjectFactory::Instance();
  factory.add(1,build1);
  factory.add(2,build2);
  factory.add(3,build3);
}

  
int main()
{
  // fill factory
  loadFactory();
  // get the factory
  auto & factory=MyObjectFactory::Instance();
  auto list =factory.registered();
  std::cout<<"Registered identifiers:"<<std::endl;
  for (std::size_t i=0; i<list.size();++i)
    {
      std::cout<<list[i]<<std::endl;
    }
  std::cout<<" Type identifier of desired object (1,2,3)"<<std::endl;
  Identifier i;
  std::cin>>i;
  auto object = factory.create(i);
  object->whoAmI();
}
