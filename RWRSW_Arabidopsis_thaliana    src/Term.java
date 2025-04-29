import java.util.ArrayList;

// 基因本体术语结构信息(项)、GO拓扑结构
public class Term {
    String id; // 基因本体术语ID

    String name; // 基因本体术语名称

    String namespace; // 基因本体术语名称所属分支

    boolean isObsolete; // 基因本体术语是否已经被淘汰

    ArrayList<String> partID = new ArrayList<>(); // part_of关联的基因本体术语ID

    ArrayList<String> isID = new ArrayList<>(); // is_a关联的基因本体术语ID

    public String getId() {
        return id;
    }

}
