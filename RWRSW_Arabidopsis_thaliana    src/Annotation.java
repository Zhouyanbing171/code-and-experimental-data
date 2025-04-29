import java.util.ArrayList;

// 注释信息
public class Annotation {

    String geneName; // 基因名称

    String geneID; // 标识符ID（与基因本体数据库中的基因对应）

    String goID; // 基因本体术语ID

    String evidenceCode; // 证据编码

    String nameSpace; // 基因本体术语所属分支

    ArrayList<String> synonym = new ArrayList<>(); // 基因本体术语的同义词

    public String getGeneName() {
        return geneName;
    }

    public String getGoID() {
        return goID;
    }

    public String getGeneID() {
        return geneID;
    }

}
