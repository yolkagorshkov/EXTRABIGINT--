#ifndef EXTRABIGINT
#define EXTRABIGINT
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
using namespace std;
class extrabigint{
public:

    friend ostream& operator<< (ostream &out, const extrabigint &ebi);
    friend istream& operator>> (istream &in, extrabigint &ebi);

    void ebi_print(string spliter=""){//очень плохое название
        extrabigint a=*this;
        int u=a.index_unzero_claster(a);
        if(a.claster[u]>0){
            for(int i=u;i<a.r;++i){
                if(i==u)cout<<a.claster[i]<<spliter;
                else{
                    string t="";
                    for(int j=0;j<18-to_string(a.claster[i]).length();++j)t+='0';
                    t+=to_string(a.claster[i]);
                    cout<<t<<spliter;
                }
            }
        }
        else{
            for(int i=u;i<a.r;++i){
                if(i==u)cout<<a.claster[i]<<spliter;
                else{
                    string t="";
                    for(int j=0;j<18-to_string(abs(a.claster[i])).length();++j)t+='0';
                    t+=to_string(abs(a.claster[i]));
                    cout<<t<<spliter;
                }
            }
        }
    }

    void ebiabs(){//еще одно плохое название
        for(int i=0;i<(*this).r;++i){
            (*this).claster[i]=abs((*this).claster[i]);
        }
    }

    void neg(){//хачу киндер пингви, злой EXTRABIGINT++ такой
        for(int i=0;i<(*this).r;++i){
            (*this).claster[i]=-abs((*this).claster[i]);
        }
    }

    extrabigint operator*(extrabigint b){
        extrabigint a=*this;
        return a.smart_multiplication(a,b);
    }

    extrabigint operator+(extrabigint b){
        extrabigint a=*this;
        return a.plus(a,b);
    }

    extrabigint operator-(extrabigint b){
        extrabigint a=*this;
        return a.minus(a,b);
    }

    extrabigint operator/(long long b){
        extrabigint a=*this;
        return a.div(b);
    }

    bool operator<(extrabigint b){
        extrabigint a=*this;
        return (a.true_compare(b)=='<');
    }

    bool operator>(extrabigint b){
        extrabigint a=*this;
        return (a.true_compare(b)=='>');
    }

    bool operator==(extrabigint b){
        extrabigint a=*this;
        return (a.true_compare(b)=='=');
    }

    bool operator!=(extrabigint b){
        extrabigint a=*this;
        return (a.true_compare(b)!='=');
    }

    long long& operator[](const int index);

    
private:

    int r=37;
    int*n=&r;
    long long*claster=new long long[*n];//самая большая моя ошибка

    int index_unzero_claster(extrabigint a){
        for(int i=0;i<a.r;i++){
            if(a.claster[i]!=0)return i;
        }
        return -1;
    }

    int strlen(extrabigint a){//надо так было для общего благополучия
        a.ebiabs();
        string as="";
        for(int i=index_unzero_claster(a);i<a.r;++i){
            string u=to_string(a.claster[i]);
            int u2=(u!="0")?u.length():0;
            if(i!=index_unzero_claster(a))for(int i=0;i<18-u2;++i)as+='0';
            as+=u;
        }
        return as.length()-a.r+index_unzero_claster(a)+1;
    }

    void init(string a){//я смотрел репортаж рен-тв во время создания этой функции
        int l=a.length();
        for(int i=l-1,sc=0,pos=(*this).r-1;i>-1;--i,++sc){
            (*this).claster[pos]+=(a[i]-'0')*pow(10,sc);
            if(sc==17){
                sc=-1;
                --pos;
            }
        }
    }

    void upper_zeroing(long long*lst,int times){//нулями все заполняется для того, чтобы не было ошибок
        for(int i=0;i<times;++i)lst[i]=0;
    }

    void del_upper_zeroing(long long*lst){// чтобы были ошибки
        for(int i=*n;i>0;--i)if(lst[i]!=000000000000000000){*n=i+1;break;}
    }

    extrabigint check(extrabigint a){//основная функция для отладки архитектуры
        long long t=pow(10,18);
        for(int i=index_unzero_claster(a);i<a.r;++i){
            if(a.claster[i]>999999999999999999){a.claster[i-1]+=a.claster[i]/t;a.claster[i]%=t;}
            if(a.claster[i]<-999999999999999999){a.claster[i-1]-=a.claster[i]/t;a.claster[i]%=t;}
        }
        for(int i=index_unzero_claster(a);i<a.r;++i){
            if(a.claster[i]>0&&a.claster[i+1]<0){
                --a.claster[i];
                a.claster[i+1]+=t;
            }
            if(a.claster[i]<0&&a.claster[i+1]>0){
                ++a.claster[i];
                a.claster[i+1]-=t;
            }
        }
        return a;
    }

    string sign_map(extrabigint a){//знаковость разрядов кластера
        string map="";
        for(int i=0;i<a.r;i++){
            if(a.claster[i]>0)map+='+';
            if(a.claster[i]<0)map+='-';
            else map+='0';
        }
        return map;
    }

    

    extrabigint plus(extrabigint a,extrabigint b){//плюс ушел в плюс, минус больше не нужен
        bool minu=0;
        if(b.claster[b.r-1]<0&&a.claster[a.r-1]<0){
            minu=1;
            a.ebiabs();b.ebiabs();
        }
        for(int i=index_unzero_claster(a);i<a.r;++i){
            a.claster[i]+=b.claster[i];
        }
        a.check(a);
        if(minu)a.neg();
        return a;
    }

    char compare(extrabigint a,extrabigint b){//Компадро камбочо
        a.ebiabs();b.ebiabs();
        if(index_unzero_claster(a)>index_unzero_claster(b))return '>';
        if(index_unzero_claster(a)<index_unzero_claster(b))return '<';
        else{
            for(int i=index_unzero_claster(a);i<a.r;i++){
                if(a.claster[i]>b.claster[i])return '>';
                if(a.claster[i]<b.claster[i])return '<';
            }
        }
        return '=';
    }

    extrabigint minus(extrabigint a,extrabigint b){//минус ушел в минус, все в плюсе
        if((a==b)){upper_zeroing(a.claster,a.r);return a;}
        b.neg();
        return a+b;
    }

    extrabigint smart_multiplication(extrabigint a,extrabigint b){//умный самый
        int cpr=0;
        if(a.claster[a.r-1]<0)cpr++;
        if(b.claster[b.r-1]<0)cpr++;
        a.ebiabs();b.ebiabs();
        if(strlen(a)>strlen(b)){
            tr:
                string multplctr="";

                for(int i=index_unzero_claster(b);i<b.r;++i){
                    string u=to_string(b.claster[i]);
                    int u2=(u!="0")?u.length():0;
                    if(i!=index_unzero_claster(b))for(int i=0;i<18-u2;++i)multplctr+='0';
                    multplctr+=u;
                }

                extrabigint pl;
                for(int j=multplctr.length()-1,c=0;j>=0;c++,--j){
                    for(int i=index_unzero_claster(a);i>=0;--i){
                        long long pre_res=a.claster[i]*(multplctr[j]-'0');
                        long long div=pow(10,18-c%18);
                        long long div2=pow(10,c%18);
                        pl.claster[i-c/18-1]+=pre_res/div;
                        pl.claster[i-c/18]+=(pre_res%div)*((div2==0)?1:div2);
                    }
                    pl.check(pl);
                }
                if(cpr==1)pl.neg();
                return pl;
            }
        else{swap(a,b); goto tr;}
        
    }

    char true_compare(extrabigint b){//истинный Компадро камбочо
        extrabigint a=*this;
        if(a.claster[a.r-1]>0&&b.claster[b.r-1]<0)return '>';
        if(a.claster[a.r-1]<0&&b.claster[b.r-1]>0)return '<';
        if(a.claster[a.r-1]>0&&b.claster[b.r-1]>0){
            if(a.r-index_unzero_claster(a)>a.r-index_unzero_claster(b))return '>';
            if(a.r-index_unzero_claster(a)<a.r-index_unzero_claster(b))return '<';
            else{
                for(int i=index_unzero_claster(a);i<a.r;i++){
                    if(a.claster[i]>b.claster[i])return '>';
                    if(a.claster[i]<b.claster[i])return '<';
                }
            }
            return '=';
        }
        else{
            a.ebiabs();b.ebiabs();
            if(a.r-index_unzero_claster(a)>b.r-index_unzero_claster(b))return '<';
            if(a.r-index_unzero_claster(a)<b.r-index_unzero_claster(b))return '>';
            else{
                for(int i=index_unzero_claster(a);i<a.r;i++){
                    if(a.claster[i]>b.claster[i])return '<';
                    if(a.claster[i]<b.claster[i])return '>';
                }
            }
            return '=';
        }
    }

    extrabigint div(long long b){//деление
        extrabigint a=*this,buff;
        int cpr=0;
        if(a.claster[0]<0)cpr++;
        if(b<0)cpr++;
        a.ebiabs();b=abs(b);
        upper_zeroing(buff.claster,37);
        if(index_unzero_claster(a)==a.r-1&&(a.claster[a.r-1]==b)){buff.claster[buff.r-1]=1; return buff;}
        if(b==1)return a;
        if(index_unzero_claster(a)==a.r-1){
            buff.claster[buff.r-1]=a.claster[a.r-1]/b;
            return buff;
        }
        for(int i=0;;++i){
            if(i%2==0){
                for(int j=index_unzero_claster(a);j<a.r;++j){
                    buff.claster[j]+=a.claster[j]/b;
                    a.claster[j]%=b;
                }
            }
            else{
                for(int j=index_unzero_claster(a);j<a.r-1;++j){
                    if(a.claster[j]>7){a.claster[j+1]+=(8*pow(10,18));a.claster[j]-=8;}
                    if(!a.claster[j])continue;
                    else{a.claster[j+1]+=(a.claster[j]*pow(10,18));a.claster[j]=0;}
                }
            }
            if(index_unzero_claster(a)==a.r-1&&a.claster[a.r-1]<b){if(cpr==1){buff.neg();}return buff;}
            if(index_unzero_claster(a)==-1){if(cpr==1){buff.neg();}return buff;}
        }
    }

};

ostream& operator<< (ostream &out, const extrabigint &a){
    
    int u=((extrabigint)a).index_unzero_claster(a);
    if(u==-1){
        out<<0;
        return out;
    }
    if(a.claster[u]>0){
        for(int i=u;i<a.r;++i){
            if(i==u)out<<a.claster[i];
            else{
                string t="";
                for(int j=0;j<18-to_string(a.claster[i]).length();++j)t+='0';
                t+=to_string(a.claster[i]);
                out<<t;
            }
        }
    }
    else{
        for(int i=u;i<a.r;++i){
            if(i==u)out<<a.claster[i];
            else{
                string t="";
                for(int j=0;j<18-to_string(abs(a.claster[i])).length();++j)t+='0';
                t+=to_string(abs(a.claster[i]));
                out<<t;
            }
        }
    }
    return out;
}

istream& operator>>(istream &in,extrabigint &ebi){
    string u;
    in>>u;
    ebi.init(u);
    return in;
}
long long& extrabigint::operator[] (const int index){
    return claster[(*this).r-1-(index%37)];
}
#endif
